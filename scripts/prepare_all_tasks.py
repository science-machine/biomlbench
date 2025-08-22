#!/usr/bin/env python3
"""
Script to prepare all BioML-bench tasks in parallel.

This script reads the production-jobs.txt file, extracts unique task IDs,
and prepares them using the biomlbench prepare command in parallel.
"""

import argparse
import concurrent.futures
import logging
import subprocess
import sys
import time
from pathlib import Path
from typing import List, Set, Tuple

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('task_preparation.log')
    ]
)
logger = logging.getLogger(__name__)


def extract_unique_tasks(jobs_file: Path) -> List[str]:
    """
    Extract unique task IDs from the production jobs file.
    
    Args:
        jobs_file: Path to the production-jobs.txt file
        
    Returns:
        List of unique task IDs
    """
    unique_tasks: Set[str] = set()
    
    try:
        with open(jobs_file, 'r') as f:
            for line in f:
                line = line.strip()
                # Skip comments and empty lines
                if not line or line.startswith('#'):
                    continue
                
                # Parse agent,task_id format
                if ',' in line:
                    _, task_id = line.split(',', 1)
                    unique_tasks.add(task_id.strip())
        
        task_list = sorted(list(unique_tasks))
        logger.info(f"Found {len(task_list)} unique tasks to prepare")
        return task_list
        
    except FileNotFoundError:
        logger.error(f"Jobs file not found: {jobs_file}")
        raise
    except Exception as e:
        logger.error(f"Error reading jobs file: {e}")
        raise


def prepare_single_task(task_id: str, extra_args: List[str] = None) -> Tuple[str, bool, str]:
    """
    Prepare a single task using biomlbench prepare command.
    
    Args:
        task_id: Task ID to prepare
        extra_args: Additional arguments to pass to biomlbench prepare
        
    Returns:
        Tuple of (task_id, success, error_message)
    """
    cmd = ['biomlbench', 'prepare', '-t', task_id]
    if extra_args:
        cmd.extend(extra_args)
    
    start_time = time.time()
    logger.info(f"Starting preparation of task: {task_id}")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=36000  # 1 hour timeout per task
        )
        
        duration = time.time() - start_time
        
        if result.returncode == 0:
            logger.info(f"✓ Successfully prepared {task_id} in {duration:.1f}s")
            return task_id, True, ""
        else:
            error_msg = f"Command failed with return code {result.returncode}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
            logger.error(f"✗ Failed to prepare {task_id} after {duration:.1f}s: {error_msg}")
            return task_id, False, error_msg
            
    except subprocess.TimeoutExpired:
        error_msg = f"Task preparation timed out after 1 hour"
        logger.error(f"✗ {task_id}: {error_msg}")
        return task_id, False, error_msg
        
    except Exception as e:
        error_msg = f"Unexpected error: {str(e)}"
        logger.error(f"✗ {task_id}: {error_msg}")
        return task_id, False, error_msg


def prepare_tasks_parallel(task_ids: List[str], num_threads: int, extra_args: List[str] = None) -> None:
    """
    Prepare multiple tasks in parallel using ThreadPoolExecutor.
    
    Args:
        task_ids: List of task IDs to prepare
        num_threads: Number of parallel threads to use
        extra_args: Additional arguments to pass to biomlbench prepare
    """
    logger.info(f"Starting parallel preparation of {len(task_ids)} tasks using {num_threads} threads")
    
    successful_tasks = []
    failed_tasks = []
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        # Submit all tasks
        future_to_task = {
            executor.submit(prepare_single_task, task_id, extra_args): task_id
            for task_id in task_ids
        }
        
        # Process completed tasks
        for future in concurrent.futures.as_completed(future_to_task):
            task_id, success, error_msg = future.result()
            
            if success:
                successful_tasks.append(task_id)
            else:
                failed_tasks.append((task_id, error_msg))
    
    # Summary
    logger.info(f"\n{'='*60}")
    logger.info(f"PREPARATION SUMMARY")
    logger.info(f"{'='*60}")
    logger.info(f"Total tasks: {len(task_ids)}")
    logger.info(f"Successful: {len(successful_tasks)}")
    logger.info(f"Failed: {len(failed_tasks)}")
    
    if successful_tasks:
        logger.info(f"\n✓ Successfully prepared tasks:")
        for task in successful_tasks:
            logger.info(f"  - {task}")
    
    if failed_tasks:
        logger.error(f"\n✗ Failed tasks:")
        for task, error in failed_tasks:
            logger.error(f"  - {task}: {error}")
        
        # Write failed tasks to file for retry
        failed_file = Path("failed_tasks.txt")
        with open(failed_file, 'w') as f:
            for task, _ in failed_tasks:
                f.write(f"{task}\n")
        logger.info(f"Failed task IDs written to: {failed_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Prepare all BioML-bench tasks in parallel",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Prepare all tasks with 4 threads
  python scripts/prepare_all_tasks.py --threads 4
  
  # Prepare with custom jobs file and keep raw data
  python scripts/prepare_all_tasks.py --jobs custom-jobs.txt --threads 8 --keep-raw
  
  # Prepare only specific tasks from a file
  python scripts/prepare_all_tasks.py --task-list failed_tasks.txt --threads 2
        """
    )
    
    parser.add_argument(
        '--jobs',
        type=Path,
        default=Path('scripts/gcp-deploy/production-jobs.txt'),
        help='Path to jobs file (default: scripts/gcp-deploy/production-jobs.txt)'
    )
    
    parser.add_argument(
        '--task-list',
        type=Path,
        help='Path to file with specific task IDs to prepare (one per line)'
    )
    
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of parallel threads (default: 4)'
    )
    
    parser.add_argument(
        '--keep-raw',
        action='store_true',
        help='Keep raw downloaded data'
    )
    
    parser.add_argument(
        '--overwrite-checksums',
        action='store_true',
        help='Overwrite existing checksums'
    )
    
    parser.add_argument(
        '--overwrite-leaderboard',
        action='store_true',
        help='Overwrite existing leaderboard'
    )
    
    parser.add_argument(
        '--skip-verification',
        action='store_true',
        help='Skip checksum verification'
    )
    
    parser.add_argument(
        '--data-dir',
        type=Path,
        help='Custom data directory'
    )
    
    args = parser.parse_args()
    
    # Validate thread count
    if args.threads < 1:
        logger.error("Number of threads must be at least 1")
        sys.exit(1)
    
    # Get task IDs
    if args.task_list:
        logger.info(f"Reading task IDs from: {args.task_list}")
        try:
            with open(args.task_list, 'r') as f:
                task_ids = [line.strip() for line in f if line.strip() and not line.startswith('#')]
        except FileNotFoundError:
            logger.error(f"Task list file not found: {args.task_list}")
            sys.exit(1)
    else:
        logger.info(f"Extracting task IDs from: {args.jobs}")
        task_ids = extract_unique_tasks(args.jobs)
    
    if not task_ids:
        logger.error("No tasks found to prepare")
        sys.exit(1)
    
    # Build extra arguments for biomlbench prepare
    extra_args = []
    if args.keep_raw:
        extra_args.append('--keep-raw')
    if args.overwrite_checksums:
        extra_args.append('--overwrite-checksums')
    if args.overwrite_leaderboard:
        extra_args.append('--overwrite-leaderboard')
    if args.skip_verification:
        extra_args.append('--skip-verification')
    if args.data_dir:
        extra_args.extend(['--data-dir', str(args.data_dir)])
    
    # Start preparation
    logger.info(f"Preparing {len(task_ids)} tasks with {args.threads} threads")
    if extra_args:
        logger.info(f"Extra arguments: {' '.join(extra_args)}")
    
    try:
        prepare_tasks_parallel(task_ids, args.threads, extra_args)
    except KeyboardInterrupt:
        logger.info("Preparation interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error during preparation: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main() 