"""
Agent execution module for BioML-bench.

This module provides functionality to run AI agents on biomedical tasks
within the biomlbench framework.
"""

import asyncio
import json
import logging
import os

# Import from agents directory
import sys
import tempfile
import time
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any, List, Tuple

import docker

from biomlbench.data import is_dataset_prepared
from biomlbench.registry import Task, registry
from biomlbench.utils import (
    create_run_dir,
    generate_submission_from_metadata,
    get_logger,
    get_runs_dir,
    get_timestamp,
)

sys.path.append(str(Path(__file__).parent.parent / "agents"))
from registry import Agent
from registry import registry as agent_registry
from run import run_in_container

from environment.defaults import DEFAULT_CONTAINER_CONFIG_PATH

logger = get_logger(__name__)


async def check_agent_success(run_dir: Path, run_logger: logging.Logger) -> bool:
    """
    Check if an agent run actually succeeded by examining submission files and logs.

    This catches cases where agents fail internally but exit with success codes.
    """
    try:
        # Check 1: Look for submission files
        submission_dir = run_dir / "submission"
        has_submission = False

        if submission_dir.exists():
            # Check for common submission file formats
            submission_formats = ["submission.csv", "submission.h5ad"]
            for format_name in submission_formats:
                submission_file = submission_dir / format_name
                if submission_file.exists() and submission_file.stat().st_size > 0:
                    has_submission = True
                    break

        # Check 2: Examine logs for failure patterns
        log_file = run_dir / "run.log"
        has_critical_errors = False

        if log_file.exists():
            try:
                with open(log_file, "r") as f:
                    log_content = f.read().lower()

                    # Critical error patterns that indicate failure
                    error_patterns = [
                        "traceback (most recent call last)",
                        "authenticationerror",
                        "invalid api key",
                    ]

                    for pattern in error_patterns:
                        if pattern in log_content:
                            has_critical_errors = True
                            run_logger.debug(f"Found error pattern '{pattern}' in logs")
                            break
            except Exception:
                run_logger.warning("Could not read run log for failure analysis")

        # # Decision logic: If we have critical errors, it's a failure
        # if has_critical_errors:
        #     run_logger.error(f"Agent failed: Critical errors found in logs")
        #     return False

        # If we have a valid submission file, it's likely a success
        if has_submission:
            run_logger.info(f"Agent succeeded: Valid submission file found")
            return True

        # If no submission and no clear errors, treat as suspicious failure
        run_logger.warning(f"Agent success unclear: No submission file found")
        return False

    except Exception as e:
        run_logger.error(f"Error during success check: {e}")
        return False


@dataclass(frozen=True)
class AgentTask:
    """Represents a single agent-task execution."""

    run_id: str
    seed: int
    image: str
    path_to_run_group: Path
    path_to_run: Path
    agent: Agent
    task: Task
    container_config: dict[str, Any]


async def worker(
    idx: int,
    queue: asyncio.Queue[AgentTask],
    client: docker.DockerClient,
    tasks_outputs: dict[str, dict[str, Any]],
    retain_container: bool = False,
) -> None:
    """Worker function that processes agent tasks from the queue."""
    while True:
        agent_task = await queue.get()

        # Create logger for the run
        run_logger = get_logger(str(agent_task.path_to_run))
        log_file_handler = logging.FileHandler(agent_task.path_to_run / "run.log")
        log_file_handler.setFormatter(logging.getLogger().handlers[0].formatter)
        run_logger.addHandler(log_file_handler)
        run_logger.propagate = False

        run_logger.info(
            f"[Worker {idx}] Running seed {agent_task.seed} for {agent_task.task.id} and agent {agent_task.agent.name}"
        )

        task_output = {
            "task_id": agent_task.task.id,
            "agent_id": agent_task.agent.id,
            "seed": agent_task.seed,
        }
        try:
            await asyncio.to_thread(
                run_in_container,
                client=client,
                task_id=agent_task.task.id,
                public_dir=agent_task.task.public_dir,
                private_dir=agent_task.task.private_dir,
                agent=agent_task.agent,
                image=agent_task.agent.name,
                container_config=agent_task.container_config,
                retain_container=retain_container,
                run_dir=agent_task.path_to_run,
                logger=run_logger,
            )

            # ENHANCED: Check if the agent actually succeeded beyond just not throwing an exception
            success = await check_agent_success(agent_task.path_to_run, run_logger)
            task_output["success"] = success

            if success:
                run_logger.info(
                    f"[Worker {idx}] Finished running seed {agent_task.seed} for {agent_task.task.id} and agent {agent_task.agent.name}"
                )
            else:
                run_logger.error(
                    f"[Worker {idx}] Agent completed without exception but failed internally for seed {agent_task.seed}, task {agent_task.task.id}, agent {agent_task.agent.name}"
                )
        except Exception as e:
            stack_trace = traceback.format_exc()
            run_logger.error(type(e))
            run_logger.error(stack_trace)
            run_logger.error(
                f"Run failed for seed {agent_task.seed}, agent {agent_task.agent.id} and task "
                f"{agent_task.task.id}"
            )
            task_output["success"] = False
        finally:
            tasks_outputs[agent_task.run_id] = task_output
            queue.task_done()


def create_task_list_file(task_id: str) -> str:
    """Create a temporary task list file for single task execution."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write(task_id + "\n")
        return f.name


def get_task_ids(task_id: str = None, task_list: str = None) -> List[str]:
    """Get list of task IDs from either single task or task list file."""
    if task_id:
        return [task_id]
    elif task_list:
        with open(task_list, "r") as f:
            return [line.strip() for line in f.read().splitlines() if line.strip()]
    else:
        raise ValueError("Either task_id or task_list must be provided")


async def run_agent_async(
    agent_id: str,
    task_ids: List[str],
    n_workers: int = 1,
    n_seeds: int = 1,
    container_config_path: str = None,
    retain_container: bool = False,
    data_dir: str = None,
) -> Tuple[str, Path]:
    """
    Run an agent on multiple tasks asynchronously.

    Returns:
        Tuple[str, Path]: The run group ID and path to the generated submission file
    """
    client = docker.from_env()

    # Set up registry with data directory
    if data_dir:
        task_registry = registry.set_data_dir(Path(data_dir))
    else:
        task_registry = registry

    # Get agent
    agent = agent_registry.get_agent(agent_id)

    # Check for privileged container requirements
    if agent.privileged and not (
        os.environ.get("I_ACCEPT_RUNNING_PRIVILEGED_CONTAINERS", "False").lower()
        in ("true", "1", "t")
    ):
        raise ValueError(
            "Agent requires running in a privileged container, but the environment variable "
            "`I_ACCEPT_RUNNING_PRIVILEGED_CONTAINERS` is not set to `True`! "
            "Carefully consider if you wish to run this agent before continuing. "
            "See agents/README.md for more details."
        )

    # Create run group
    run_group = f"{get_timestamp()}_run-group_{agent.name}"

    # Validate all tasks are prepared
    for task_id in task_ids:
        task = task_registry.get_task(task_id)
        if not is_dataset_prepared(task):
            raise ValueError(
                f"Dataset for task `{task.id}` is not prepared! "
                f"Please run `biomlbench prepare -t {task.id}` to prepare the dataset."
            )

    # Load container configuration
    if container_config_path is None:
        container_config_path = DEFAULT_CONTAINER_CONFIG_PATH

    with open(container_config_path, "r") as f:
        container_config = json.load(f)

    # Create agent tasks for each (task × seed) combination
    logger.info(f"Launching run group: {run_group}")
    agent_tasks = []
    for seed in range(n_seeds):
        for task_id in task_ids:
            task = task_registry.get_task(task_id)
            run_dir = create_run_dir(task.id, agent.id, run_group)
            # Store relative path from run group directory for proper reconstruction
            run_group_dir = get_runs_dir() / run_group
            run_id = str(run_dir.relative_to(run_group_dir))
            agent_task = AgentTask(
                run_id=run_id,
                seed=seed,
                image=agent.name,
                agent=agent,
                task=task,
                path_to_run_group=run_dir.parent,
                path_to_run=run_dir,
                container_config=container_config,
            )
            agent_tasks.append(agent_task)

    logger.info(f"Creating {n_workers} workers to serve {len(agent_tasks)} tasks...")

    # Create queue and workers
    queue = asyncio.Queue()
    for agent_task in agent_tasks:
        queue.put_nowait(agent_task)

    workers = []
    tasks_outputs = {}
    for idx in range(n_workers):
        w = asyncio.create_task(worker(idx, queue, client, tasks_outputs, retain_container))
        workers.append(w)

    # Wait for completion
    started_at = time.monotonic()
    await queue.join()
    time_taken = time.monotonic() - started_at

    # Clean up workers
    for w in workers:
        w.cancel()
    await asyncio.gather(*workers, return_exceptions=True)

    # Generate metadata
    metadata = {
        "run_group": run_group,
        "created_at": get_timestamp(),
        "runs": tasks_outputs,
        "agent_id": agent_id,
        "task_ids": task_ids,
        "n_seeds": n_seeds,
        "n_workers": n_workers,
    }

    run_group_dir = get_runs_dir() / run_group
    metadata_path = run_group_dir / "metadata.json"
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=4, sort_keys=False, default=str)

    # Auto-generate submission file for grading
    submission_path = generate_submission_from_metadata(metadata_path)

    # CRITICAL: Check for failures and report them to prevent silent failures
    failed_runs = [
        run_id for run_id, output in tasks_outputs.items() if not output.get("success", False)
    ]
    total_runs = len(tasks_outputs)
    successful_runs = total_runs - len(failed_runs)

    logger.info(f"{n_workers} workers ran for {time_taken:.2f} seconds in total")
    logger.info(f"Results saved to: {run_group_dir}")

    # Report success/failure statistics
    if failed_runs:
        logger.error(f"❌ {len(failed_runs)}/{total_runs} runs FAILED!")
        logger.error(f"✅ {successful_runs}/{total_runs} runs succeeded")

        # Check for common failure patterns and provide specific guidance
        first_failed_run = failed_runs[0]
        first_run_log = run_group_dir / first_failed_run / "run.log"

        if first_run_log.exists():
            try:
                with open(first_run_log, "r") as f:
                    log_content = f.read().lower()

                    if (
                        "docker.errors.imagenotfound" in log_content
                        or "no such image" in log_content
                    ):
                        logger.error(f"🐳 Docker image '{agent.name}:latest' not found!")
                        logger.error(f"💡 To fix this, build the agent image first:")
                        logger.error(f"   ./scripts/build_agent.sh {agent.name}")
                        logger.error(f"   # OR build all agent images:")
                        logger.error(f"   ./scripts/build_agent.sh --all")
                    elif "environment variable" in log_content and "not set" in log_content:
                        logger.error(f"🔑 Missing API key or environment variable!")
                        logger.error(f"💡 To fix this, check your .env file in the project root")
                        logger.error(f"   Example: OPENAI_API_KEY=your-key-here")
                    else:
                        logger.error(f"📋 Check detailed error logs:")
                        logger.error(f"   cat {first_run_log}")
            except Exception:
                logger.error(f"📋 Check detailed error logs in: {run_group_dir}")

        logger.error(f"🚨 AGENT EXECUTION FAILED - See errors above")

        # If ALL runs failed, this is a critical error that should not be silent
        if successful_runs == 0:
            raise RuntimeError(
                f"All {total_runs} agent runs failed! "
                f"This indicates a systematic issue (missing Docker image, API keys, etc.). "
                f"Check the error messages above and the detailed logs in: {run_group_dir}"
            )
    else:
        logger.info(f"✅ All {total_runs} runs completed successfully!")

    if successful_runs > 0:
        logger.info(f"Submission file ready for grading: {submission_path}")
        logger.info(
            f"To grade results, run: biomlbench grade --submission {submission_path} --output-dir results/"
        )
    else:
        logger.warning(
            "⚠️  No successful runs to grade - submission file contains no valid results"
        )

    return run_group, submission_path


def run_agent(args) -> str:
    """
    Main entry point for running agents from the CLI.

    Args:
        args: Parsed command line arguments

    Returns:
        str: The run group ID for this execution
    """
    # Get task IDs
    task_ids = get_task_ids(args.task_id, args.task_list)

    logger.info(f"Running agent '{args.agent}' on tasks: {task_ids}")

    # Run the agent
    run_group, submission_path = asyncio.run(
        run_agent_async(
            agent_id=args.agent,
            task_ids=task_ids,
            n_workers=args.n_workers,
            n_seeds=args.n_seeds,
            container_config_path=args.container_config,
            retain_container=args.retain_container,
            data_dir=args.data_dir,
        )
    )

    return run_group
