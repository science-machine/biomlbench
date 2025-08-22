import logging
import time
from pathlib import Path

import docker
from docker.models.containers import Container
from dotenv import dotenv_values

from agents.registry import Agent
from biomlbench.registry import Task
from biomlbench.utils import purple
from environment.utils import (
    create_task_container,
    extract_from_container,
    extract_from_container_sysbox,
)

CONSTANTS = dotenv_values(Path(__file__).parent.resolve() / ".shared_env")


def save_output(container: Container, save_dir: Path, container_config: dict) -> Path:
    """
    Extracts the submission, logs, and code directories from the container

    and saves them to the specified directory.

    Args:
        container: The Docker container.
        save_dir: The directory where the output file will be saved.
        container_config: The container configuration.
    Returns:
        Path to the output directory.
    """
    if "runtime" in container_config and container_config["runtime"] == "sysbox-runc":
        extraction_fn = extract_from_container_sysbox
    else:
        extraction_fn = extract_from_container

    for dir_type in ["SUBMISSION_DIR", "LOGS_DIR", "CODE_DIR"]:
        container_dir = CONSTANTS[dir_type]
        extraction_fn(container, container_dir, save_dir)

    return save_dir


def execute_agent(container: Container, agent: Agent, logger: logging.Logger):
    """
    Initiates the agent via its start script inside the container.
    """
    cmd = ["bash", f"{CONSTANTS['AGENT_DIR']}/start.sh"]

    if agent.kwargs_type == "argparse":
        for key, value in agent.kwargs.items():
            cmd += [f"--{key}", str(value)]

    if agent.kwargs_type == "omegaconf":
        cmd += [f"{key}={value}" for key, value in agent.kwargs.items()]

    logger.info("Running agent...")
    exit_code, output = container.exec_run(cmd, stream=True, user="nonroot")

    for chunk in output:
        logger.info(f"[Container] {chunk.decode('utf-8').strip()}")


def clean_up(container: Container, logger: logging.Logger, retain: bool = False) -> bool:
    """
    Stops and removes the container.

    Returns:
        True if successful, False otherwise.
    """
    logger.info(f"Cleaning up container: {container.name}")
    try:
        container.stop()
        if not retain:
            container.remove()
        logger.info(f"Container {container.name} stopped and removed.")
        return True
    except Exception as e:
        logger.error(
            f"Error cleaning up: {e}. You may wish to manually check the status of the {container.name} container."
        )
        return False


def run_in_container(
    client: docker.DockerClient,
    task_id: str,
    public_dir: Path,
    private_dir: Path,
    agent: Agent,
    image: str,
    container_config: dict,
    retain_container: bool,
    run_dir: Path,
    logger: logging.Logger,
) -> Path:
    """
    Runs environment containing the task and agent for a set maximum amount of time.

    Args:
        client: Docker client.
        task_id: The task to run.
        public_dir: The public directory of the task.
        private_dir: The private directory of the task.
        agent: The agent to run.
        image: The Docker image to use. Assumes the image is built.
        container_config: Configuration for the Docker container.
        retain_container: Whether to retain the container after the run instead of removing it.
        run_dir: Path to the directory where all assets associated with the run are stored.
        logger: Logger for the run.

    Returns:
        Path to the output file.
    """
    volumes_config = {
        public_dir.resolve().as_posix(): {
            "bind": "/home/data",
            "mode": "ro",
        },
        private_dir.resolve().as_posix(): {
            "bind": f"/private/data/{task_id}/prepared/private/",
            "mode": "ro",
        },
    }

    container = create_task_container(
        client=client,
        task_id=task_id,
        container_config=container_config,
        volumes_config=volumes_config,
        env_vars={
            "TASK_ID": task_id,
            **agent.env_vars,
        },
        container_image=image,
        privileged=agent.privileged,
    )

    logger.info(purple(f"Run started: {run_dir}"))
    try:
        time_start = time.monotonic()
        container.start()
        logger.info("Waiting for grading server to start...")
        
        # Add more verbose health check with actual output
        logger.info("Running health check...")
        TIMEOUT_TIME = 2400
        exit_code, output = container.exec_run(
            f'timeout {TIMEOUT_TIME}s sh -c "while ! curl -s http://localhost:5000/health > /dev/null; do sleep 1; done"'
        )
        
        if exit_code != 0:
            # Get more detailed error information
            logger.error(f"Health check failed with exit code {exit_code}")
            
            # Try to get container logs
            try:
                container_logs = container.logs(stdout=True, stderr=True, tail=50).decode('utf-8')
                logger.error("Container logs (last 50 lines):")
                for line in container_logs.split('\n'):
                    if line.strip():
                        logger.error(f"[Container] {line}")
            except Exception as log_error:
                logger.error(f"Could not retrieve container logs: {log_error}")
            
            # Try to check if entrypoint.sh is running
            try:
                ps_exit, ps_output = container.exec_run("ps aux")
                if ps_exit == 0:
                    logger.error("Running processes in container:")
                    for line in ps_output.decode('utf-8').split('\n'):
                        if line.strip():
                            logger.error(f"[Process] {line}")
                else:
                    logger.error("Could not get process list from container")
            except Exception as ps_error:
                logger.error(f"Could not get process info: {ps_error}")
            
            # Try to check if entrypoint log exists
            try:
                log_check_exit, log_check_output = container.exec_run("ls -la /home/logs/ || echo 'No logs dir'")
                logger.error("Contents of /home/logs/:")
                logger.error(log_check_output.decode('utf-8'))
                
                # Try to read entrypoint log if it exists
                entrypoint_log_exit, entrypoint_log_output = container.exec_run("tail -50 /home/logs/entrypoint.log 2>/dev/null || echo 'No entrypoint.log'")
                if entrypoint_log_exit == 0:
                    logger.error("Entrypoint log (last 50 lines):")
                    for line in entrypoint_log_output.decode('utf-8').split('\n'):
                        if line.strip():
                            logger.error(f"[Entrypoint] {line}")
            except Exception as entrypoint_error:
                logger.error(f"Could not check entrypoint log: {entrypoint_error}")
            
            raise RuntimeError(
                f"The grading server failed to start within {TIMEOUT_TIME} seconds. This is likely due to an error in `entrypoint.sh`; check the logs above for details."
            )
            
        execute_agent(container, agent, logger)
        save_output(container, run_dir, container_config)
        time_end = time.monotonic()
        logger.info(f"Run completed in {time_end - time_start:.2f} seconds.")
        return run_dir
    except Exception as e:
        raise e
    finally:
        clean_up(container, logger, retain_container)
