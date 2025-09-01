#!/usr/bin/env python3
"""
Local BioML-bench Test Script

Simple test runner that:
1. Reads jobs from a text file
2. Runs biomlbench pipeline locally for each task
3. Manages parallel execution
4. Provides progress tracking and logging
"""

import argparse
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Tuple


def log(message: str):
    """Simple logging with timestamp."""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")


def run_command(
    cmd: List[str], description: str = "", timeout: int = None, cwd: str = None
) -> Tuple[int, str]:
    """Run a command and return (exit_code, output)."""
    try:
        if description:
            log(f"Running: {description}")

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, cwd=cwd)
        return result.returncode, result.stdout + result.stderr
    except subprocess.TimeoutExpired:
        return 124, "Command timed out"
    except Exception as e:
        return 1, str(e)


def run_biomlbench_job(agent: str, task_id: str, project_root: str) -> bool:
    """Run the biomlbench pipeline locally."""
    log(f"Running job locally: {agent} -> {task_id}")

    try:
        # Activate virtual environment and run agent
        cmd = [
            "bash",
            "-c",
            f"""
            set -e
            cd {project_root}
            source .venv/bin/activate
            
            echo "📋 Running agent: {agent} on task: {task_id}"
            biomlbench run-agent --agent {agent} --task-id {task_id}
            
            # Get the specific run group ID that was just created (now includes VM UUID if provided)
            LATEST_RUN=$(find runs/ -name "*run-group_{agent}*" -type d | sort | tail -1)
            if [ -z "$LATEST_RUN" ]; then
                echo "❌ No run directory found for agent {agent}"
                exit 1
            fi
            
            RUN_GROUP_ID=$(basename "$LATEST_RUN")
            echo "📁 Run group: $RUN_GROUP_ID"
            
            # Grade results
            echo "📊 Grading results..."
            biomlbench grade --submission "$LATEST_RUN/submission.jsonl" --output-dir results/
            
            # Get the grading timestamp from the most recent grading report
            GRADING_REPORT=$(find results/ -name "*_grading_report.json" | sort | tail -1)
            if [ -n "$GRADING_REPORT" ]; then
                GRADING_TIMESTAMP=$(basename "$GRADING_REPORT" | cut -d'_' -f1)
                echo "📊 Grading timestamp: $GRADING_TIMESTAMP"
                echo "📄 Grading report: $GRADING_REPORT"
            else
                echo "⚠️  No grading report found"
            fi
            
            echo "✅ Job completed successfully"
            """,
        ]

        exit_code, output = run_command(cmd, f"Run {agent} -> {task_id}")

        if exit_code == 0:
            log(f"✅ Job completed: {agent} -> {task_id}")
            # Show abbreviated output
            if output.strip():
                # Show only the important parts of the output
                lines = output.strip().split("\n")
                important_lines = [
                    line
                    for line in lines
                    if any(marker in line for marker in ["📋", "📁", "📊", "📄", "✅", "❌", "⚠️"])
                ]
                if important_lines:
                    print("   " + "\n   ".join(important_lines))
            return True
        else:
            log(f"❌ Job failed: {agent} -> {task_id}")
            # Show error output
            if output.strip():
                print(f"   Error output: {output.strip()}")
            return False

    except Exception as e:
        log(f"❌ Exception running job {agent} -> {task_id}: {e}")
        return False


def process_job(job: Tuple[str, str], project_root: str) -> Tuple[bool, str, str]:
    """Process a single job (agent, task_id) locally."""
    agent, task_id = job

    try:
        success = run_biomlbench_job(agent, task_id, project_root)
        return success, agent, task_id

    except Exception as e:
        log(f"❌ Exception processing job {agent} -> {task_id}: {e}")
        return False, agent, task_id


def load_jobs(jobs_file: str) -> List[Tuple[str, str]]:
    """Load jobs from a text file. Handles paths relative to current working directory."""
    jobs = []

    # Handle both absolute and relative paths
    jobs_path = Path(jobs_file)
    if not jobs_path.is_absolute():
        # If relative, try current working directory first, then relative to this script
        if (Path.cwd() / jobs_path).exists():
            jobs_path = Path.cwd() / jobs_path
        else:
            jobs_path = Path(__file__).parent / jobs_path

    try:
        with open(jobs_path, "r") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                parts = line.split(",")
                if len(parts) != 2:
                    log(f"⚠️  Skipping invalid job on line {line_num}: {line}")
                    continue

                agent, task_id = parts[0].strip(), parts[1].strip()
                jobs.append((agent, task_id))

        log(f"Loaded {len(jobs)} jobs from {jobs_path}")
        return jobs
    except FileNotFoundError:
        log(f"❌ Jobs file not found: {jobs_file}")
        log(f"   Looked for: {jobs_path}")
        return []


def find_project_root() -> str:
    """Find the biomlbench project root directory."""
    current = Path.cwd()

    # Look for common project indicators
    while current != current.parent:
        if (current / "biomlbench").exists() or (current / "pyproject.toml").exists():
            return str(current)
        current = current.parent

    # Default to current directory
    return str(Path.cwd())


def main():
    parser = argparse.ArgumentParser(
        description="Run BioML-bench jobs locally for testing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # From project root:
  python scripts/gcp-deploy/local-test.py --jobs scripts/gcp-deploy/gcp-test3.txt --concurrent 4
  
  # From scripts/gcp-deploy directory:
  python local-test.py --jobs gcp-test3.txt --concurrent 4
  
Jobs file format (one per line):
  agent,task_id
  dummy,polarishub/tdcommons-caco2-wang
  dummy,proteingym-dms/A0A1I9GEU1_NEIME
        """,
    )

    parser.add_argument("--jobs", required=True, help="Path to jobs file (agent,task_id per line)")
    parser.add_argument(
        "--concurrent", type=int, default=4, help="Maximum concurrent jobs (default: 4)"
    )
    parser.add_argument(
        "--project-root", help="Path to biomlbench project root (auto-detected if not specified)"
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="Show jobs that would be run without executing"
    )
    parser.add_argument(
        "--sequential", action="store_true", help="Run jobs sequentially instead of in parallel"
    )

    args = parser.parse_args()

    # Find project root
    project_root = args.project_root or find_project_root()
    if not os.path.exists(project_root):
        log(f"❌ Project root not found: {project_root}")
        return 1

    log(f"Using project root: {project_root}")

    # Check if virtual environment exists
    venv_path = Path(project_root) / ".venv" / "bin" / "activate"
    if not venv_path.exists():
        log(f"❌ Virtual environment not found at {venv_path}")
        log("   Please ensure .venv exists in the project root")
        return 1

    # Load jobs
    jobs = load_jobs(args.jobs)
    if not jobs:
        log("❌ No valid jobs found")
        return 1

    # Dry run
    if args.dry_run:
        log("Dry run - would execute these jobs:")
        for i, (agent, task_id) in enumerate(jobs, 1):
            print(f"  {i:2d}. {agent} -> {task_id}")
        return 0

    # Run jobs
    concurrent_jobs = 1 if args.sequential else args.concurrent
    log(
        f"Starting local test run with {concurrent_jobs} concurrent job{'s' if concurrent_jobs > 1 else ''}"
    )
    log(f"Jobs to process: {len(jobs)}")

    successful_jobs = 0
    failed_jobs = 0
    start_time = time.time()

    if args.sequential:
        # Sequential execution
        for i, job in enumerate(jobs, 1):
            agent, task_id = job
            log(f"Processing job {i}/{len(jobs)}: {agent} -> {task_id}")

            success, _, _ = process_job(job, project_root)
            if success:
                successful_jobs += 1
                log(f"🎉 SUCCESS ({i}/{len(jobs)}): {agent} -> {task_id}")
            else:
                failed_jobs += 1
                log(f"💥 FAILED ({i}/{len(jobs)}): {agent} -> {task_id}")
    else:
        # Parallel execution
        with ThreadPoolExecutor(max_workers=args.concurrent) as executor:
            # Submit all jobs
            future_to_job = {executor.submit(process_job, job, project_root): job for job in jobs}

            # Process completed jobs
            completed = 0
            for future in as_completed(future_to_job):
                completed += 1
                job = future_to_job[future]
                agent, task_id = job

                try:
                    success, _, _ = future.result()
                    if success:
                        successful_jobs += 1
                        log(f"🎉 SUCCESS ({completed}/{len(jobs)}): {agent} -> {task_id}")
                    else:
                        failed_jobs += 1
                        log(f"💥 FAILED ({completed}/{len(jobs)}): {agent} -> {task_id}")
                except Exception as e:
                    failed_jobs += 1
                    log(f"💥 EXCEPTION ({completed}/{len(jobs)}): {agent} -> {task_id}: {e}")

    # Summary
    elapsed_time = time.time() - start_time
    total_jobs = successful_jobs + failed_jobs
    log("=" * 60)
    log("LOCAL TEST COMPLETE")
    log(f"Total jobs: {total_jobs}")
    log(f"Successful: {successful_jobs}")
    log(f"Failed: {failed_jobs}")
    log(f"Success rate: {successful_jobs/total_jobs*100:.1f}%" if total_jobs > 0 else "N/A")
    log(f"Elapsed time: {elapsed_time:.1f} seconds")
    log(f"Average time per job: {elapsed_time/total_jobs:.1f} seconds" if total_jobs > 0 else "N/A")

    if successful_jobs > 0:
        log(f"Results saved in: {project_root}/runs/")
        log(f"Grading reports in: {project_root}/results/")

    return 0 if failed_jobs == 0 else 1


if __name__ == "__main__":
    exit(main())
