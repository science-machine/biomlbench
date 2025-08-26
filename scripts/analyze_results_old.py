#!/usr/bin/env python3
"""
BioML-bench Results Analysis Script

Downloads and analyzes biomlbench results from the organized S3 structure.

Features:
- Parallel downloads (configurable workers)
- Filters for organized structure only (agent/task/file)
- Agent filtering (--agents stella mlagentbench)
- Task filtering (--tasks polarishub kaggle)
- Excludes failed runs due to API quota/credit issues
- Fixes submission paths and regrades label-projection tasks
- Comprehensive analysis outputs:
  * Markdown summary report
  * CSV files for further analysis
  * Agent performance comparison
  * Task difficulty rankings
  * Detailed results matrix

Usage Examples:
  # Full analysis
  python scripts/analyze_results.py --output-dir results_analysis
  
  # Filter specific agents
  python scripts/analyze_results.py --agents stella biomni --parallel-downloads 15
  
  # Filter specific task types
  python scripts/analyze_results.py --tasks polarishub manual --output-dir polaris_manual_analysis
  
  # Dry run to see what would be downloaded
  python scripts/analyze_results.py --dry-run
"""

import argparse
import gzip
import json
import shutil
import subprocess
import tarfile
# Removed defaultdict - explicit error handling only
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Tuple, Any
import pandas as pd
import os
from scipy.stats import rankdata

def log(message: str):
    """Simple logging with timestamp."""
    import time
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def run_command(cmd: List[str], description: str = "") -> Tuple[int, str]:
    """Run a command and return (exit_code, output)."""
    if description:
        log(f"Running: {description}")
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode, result.stdout + result.stderr

def check_for_api_quota_errors(run_log_path: Path) -> bool:
    """
    Check if a run.log file contains API quota/credit error messages.
    
    Returns True if the run should be excluded (contains quota/credit errors).
    """
    if not run_log_path.exists():
        return False
    
    with open(run_log_path, 'r', encoding='utf-8', errors='ignore') as f:
        log_content = f.read()
    
    # Check for specific error phrases
    quota_phrases = [
        "credit balance",
        "quota"
    ]
    
    for phrase in quota_phrases:
        if phrase in log_content:
            log(f"🚫 Excluding failed run due to API quota/credit issue: {run_log_path}")
            return True
    
    return False

def check_for_quick_failure(run_log_path: Path, agent_id: str) -> bool:
    """
    Check if a run failed too quickly (< 600 seconds), indicating a system issue.
    
    Returns True if the run should be excluded (failed too quickly).
    Exception: dummy agent runs are never excluded for being too quick.
    """
    if agent_id == 'dummy':
        return False  # Never exclude dummy agent runs for being too quick
    
    if not run_log_path.exists():
        return False
    
    try:
        with open(run_log_path, 'r', encoding='utf-8', errors='ignore') as f:
            log_content = f.read()
    except Exception as e:
        log(f"⚠️  Error reading run.log: {e}")
        return False
    
    # Look for the completion time line - handle multiple possible formats
    import re
    completion_patterns = [
        r'\[run\.py:\d+\]\s*Run completed in ([\d.]+) seconds',  # Format with [run.py:204]
        r'Run completed in ([\d.]+) seconds',  # Simple format
        r'completed in ([\d.]+) seconds'  # Even simpler format
    ]
    
    for pattern in completion_patterns:
        match = re.search(pattern, log_content, re.IGNORECASE)
        if match:
            runtime_seconds = float(match.group(1))
            if runtime_seconds < 300:
                log(f"🚫 Excluding run due to quick failure ({runtime_seconds:.1f}s < 300s): {run_log_path}")
                return True
            break
    
    return False

def check_for_max_steps_error(run_log_path: Path, agent_id: str) -> bool:
    """
    Check if a STELLA run failed due to reaching max steps.
    
    Returns True if the run should be excluded (reached max steps).
    Only applies to STELLA agent.
    """
    if agent_id != 'stella':
        return False  # Only check STELLA runs
    
    if not run_log_path.exists():
        return False
    
    try:
        with open(run_log_path, 'r', encoding='utf-8', errors='ignore') as f:
            log_content = f.read()
    except Exception as e:
        log(f"⚠️  Error reading run.log: {e}")
        return False
    
    # Check for "Reached max steps" error
    if "Reached max steps" in log_content:
        log(f"🚫 Excluding STELLA run due to max steps error: {run_log_path}")
        return True
    
    return False

def check_for_grading_server_failure(run_log_path: Path) -> bool:
    """
    Check if a run failed due to grading server failure.
    
    Returns True if the run should be excluded (grading server failed to start).
    """
    if not run_log_path.exists():
        return False
    
    try:
        with open(run_log_path, 'r', encoding='utf-8', errors='ignore') as f:
            log_content = f.read()
    except Exception as e:
        log(f"⚠️  Error reading run.log: {e}")
        return False
    
    # Check for "The grading server failed to start" error
    if "The grading server failed to start" in log_content:
        log(f"🚫 Excluding run due to grading server failure: {run_log_path}")
        return True
    
    return False

def fix_submission_paths_and_regrade(runs_dir: Path, output_dir: Path):
    """Fix submission paths in submission.jsonl files for label-projection tasks and regrade them."""
    log("🔧 Fixing submission paths for label-projection tasks...")
    
    # Find all submission.jsonl files for label-projection tasks
    label_projection_submissions = []
    for submission_file in runs_dir.rglob('submission.jsonl'):
        with open(submission_file, 'r') as f:
            submission_data = json.load(f)
        
        if submission_data.get('task_id') == 'manual/open-problems-label-projection':
            label_projection_submissions.append((submission_file, submission_data))
    
    log(f"Found {len(label_projection_submissions)} label-projection submissions to fix")
    
    if not label_projection_submissions:
        return
    
    fixed_count = 0
    regrade_paths = []
    
    for submission_file, submission_data in label_projection_submissions:
        # Parse the original paths
        original_submission_path = submission_data['submission_path']
        original_logs_path = submission_data['logs_path']
        original_code_path = submission_data['code_path']
        
        # Extract the run-specific directory structure from the original path
        # Original: /home/runner/biomlbench/runs/2025-08-23T04-48-11-GMT_run-group_biomni/manual/open-problems-label-projection_c2555573-0104-466f-9d51-75eb4238b549/submission/submission.csv
        # New: /home/paperspace/biomlbench/analysis_resultsv2/runs/biomni/manual-open-problems-label-projection/2025-08-23T04-48-11-GMT_run-group_biomni/manual/open-problems-label-projection_c2555573-0104-466f-9d51-75eb4238b549/submission/submission.csv
        
        # Get the local base directory where this submission.jsonl exists
        local_base = submission_file.parent.resolve()
        
        # Find the actual submission.csv file locally
        actual_submission_path = None
        actual_logs_path = None
        actual_code_path = None
        
        # Look for submission/submission.csv in the local directory structure
        for potential_sub in local_base.rglob('submission.csv'):
            if potential_sub.parent.name == 'submission':
                actual_submission_path = potential_sub.resolve()
                # Get the run directory (parent of submission directory)
                run_dir = potential_sub.parent.parent.resolve()
                actual_logs_path = run_dir / 'logs'
                actual_code_path = run_dir / 'code'
                break
        
        if not actual_submission_path or not actual_submission_path.exists():
            log(f"⚠️  Could not find actual submission.csv for {submission_file}")
            continue
        
        # Update the submission data with corrected absolute paths
        submission_data['submission_path'] = str(actual_submission_path.resolve())
        submission_data['logs_path'] = str(actual_logs_path.resolve())
        submission_data['code_path'] = str(actual_code_path.resolve())
        
        # Write back the corrected submission.jsonl
        with open(submission_file, 'w') as f:
            json.dump(submission_data, f)
        
        log(f"✅ Fixed paths in {submission_file}")
        fixed_count += 1
        regrade_paths.append(submission_file)
    
    log(f"Fixed {fixed_count} submission.jsonl files")
    
    # Now regrade all the fixed submissions
    if regrade_paths:
        log(f"🎯 Regrading {len(regrade_paths)} label-projection submissions...")
        
        for submission_file in regrade_paths:
            log(f"Regrading {submission_file}")
            
            # Create output directory for this grading
            grade_output_dir = output_dir / 'regraded_label_projection' / submission_file.parent.name
            grade_output_dir.mkdir(parents=True, exist_ok=True)
            
            # Run biomlbench grade command
            cmd = [
                'uv', 'run', 'biomlbench', 'grade',
                '--submission', str(submission_file),
                '--output-dir', str(grade_output_dir)
            ]
            print(cmd)
            
            exit_code, output = run_command(cmd, f"Regrade {submission_file.name}")
            
            if exit_code == 0:
                log(f"✅ Successfully regraded {submission_file.name}")
            else:
                log(f"❌ Failed to regrade {submission_file.name}: {output}")
        
        log("🎉 Regrading complete!")

def discover_s3_artifacts_single_version(bucket: str, prefix: str, agents=None, tasks=None, exclude_tasks=None) -> Dict[str, List[str]]:
    """Discover S3 artifacts from a single version with optional task exclusion."""
    artifacts = {
        'runs': [],
        'grades': [],
        'failed_runs': [],
        'failed_grades': []
    }
    
    for artifact_type in artifacts.keys():
        s3_path = f"s3://{bucket}/{prefix}/{artifact_type}/"
        exit_code, output = run_command([
            "aws", "s3", "ls", s3_path, "--recursive"
        ], f"List {artifact_type} from {prefix}")
        
        if exit_code == 0:
            for line in output.strip().split('\n'):
                if line.strip():
                    # Parse S3 ls output: "2025-08-20 01:10:23    123 path/to/file"
                    parts = line.split()
                    if len(parts) >= 4:
                        s3_key = ' '.join(parts[3:])
                        full_path = f"s3://{bucket}/{s3_key}"
                        
                        # Check if this is organized structure or flat structure
                        key_parts = s3_key.split('/')
                        should_include = False
                        
                        if len(key_parts) >= 5:  # Potential organized structure: prefix/artifact_type/agent/task/file
                            agent_part = key_parts[-3]
                            task_part = key_parts[-2] 
                            filename = key_parts[-1]
                            
                            # Check for organized structure (agent/task/file)
                            if (not agent_part.startswith('2025-') and 
                                not task_part.startswith('2025-') and
                                filename.startswith('2025-')):
                                
                                # Apply agent filter
                                if agents and agent_part not in agents:
                                    continue
                                    
                                # Apply task filter
                                if tasks:
                                    task_matches = any(pattern in task_part for pattern in tasks)
                                    if not task_matches:
                                        continue
                                
                                # Apply exclude_tasks filter
                                if exclude_tasks:
                                    task_excluded = any(pattern in task_part for pattern in exclude_tasks)
                                    if task_excluded:
                                        continue
                                
                                should_include = True
                        
                        # For runs: include organized structure
                        if artifact_type == 'runs':
                            if should_include:
                                artifacts[artifact_type].append(full_path)
                        
                        # For grades: only include individual reports (skip aggregated reports)
                        elif artifact_type == 'grades':
                            if should_include and 'individual_reports' in filename:
                                artifacts[artifact_type].append(full_path)
                        
                        # For failed_runs/failed_grades: include both organized and flat structure
                        elif artifact_type in ['failed_runs', 'failed_grades']:
                            if should_include:
                                artifacts[artifact_type].append(full_path)
                            else:
                                # Check for flat structure: artifacts/failed_runs/2025-08-19T00-19-11-GMT_run-group_dummy.tar.gz
                                filename = key_parts[-1]
                                if filename.startswith('2025-') and '_run-group_' in filename:
                                    # Apply agent filter for flat structure
                                    if agents and '_run-group_' in filename:
                                        file_agent = filename.split('_run-group_')[-1].replace('.tar.gz', '').replace('.json.gz', '')
                                        if file_agent not in agents:
                                            continue
                                    
                                    # For flat structure failed runs, we can't easily determine the task type
                                    # So we need to be more careful about excluding manual tasks
                                    # Since we can't reliably filter flat structure by task, we'll include them
                                    # but the count_total_run_attempts function should handle the filtering
                                    artifacts[artifact_type].append(full_path)
    
    return artifacts


def discover_s3_artifacts(bucket: str='biomlbench', prefix: str='v1', agents = None, tasks= None) -> Dict[str, List[str]]:
    """Discover S3 artifacts with smart filtering, combining v1 and v2.
    
    Uses v2 artifacts for single cell (manual) tasks and v1 for everything else.
    """
    log("Discovering S3 artifacts (v1 for non-single-cell, v2 for single-cell)...")
    
    # Get v1 artifacts excluding single cell (manual) tasks
    log("Fetching v1 artifacts (excluding single cell tasks)...")
    
    # Determine what tasks to exclude from v1
    # If no specific tasks requested, exclude manual
    # If specific tasks requested and they include manual, we still need to exclude manual from v1
    v1_exclude_tasks = ['manual']
    v1_tasks = tasks
    if tasks and any('manual' in t for t in tasks):
        # If manual tasks are requested, we need to get them from v2 only
        # So remove manual from v1 tasks
        v1_tasks = [t for t in tasks if 'manual' not in t] if tasks else None
        
    v1_artifacts = discover_s3_artifacts_single_version(
        bucket=bucket, 
        prefix='v1/artifacts',
        agents=agents,
        tasks=v1_tasks,
        exclude_tasks=v1_exclude_tasks  # Always exclude manual tasks from v1
    )
    
    # Get v2 artifacts for single cell (manual) tasks only
    log("Fetching v2 artifacts (single cell tasks only)...")
    v2_artifacts = discover_s3_artifacts_single_version(
        bucket=bucket,
        prefix='v2/artifacts', 
        agents=agents,
        tasks=['manual'] if not tasks or any('manual' in t for t in tasks) else []  # Only include manual tasks from v2
    )
    
    # Combine the artifacts
    combined_artifacts = {
        'runs': [],
        'grades': [],
        'failed_runs': [],
        'failed_grades': []
    }
    
    # For v1, we already excluded 'manual' tasks, so we can add all v1 artifacts
    # For v2, we only included 'manual' tasks, so we can add all v2 artifacts
    # This ensures no overlap and that single cell (manual) tasks use ONLY v2 data
    for artifact_type in combined_artifacts.keys():
        combined_artifacts[artifact_type] = v1_artifacts[artifact_type] + v2_artifacts[artifact_type]
    
    log(f"Combined artifacts:")
    log(f"  - v1 (non-single-cell): {len(v1_artifacts['runs'])} runs, {len(v1_artifacts['grades'])} grades")
    log(f"  - v2 (single-cell): {len(v2_artifacts['runs'])} runs, {len(v2_artifacts['grades'])} grades")
    log(f"  - Total: {len(combined_artifacts['runs'])} organized run artifacts, {len(combined_artifacts['grades'])} organized grade artifacts")
    log(f"  - Total: {len(combined_artifacts['failed_runs'])} failed runs, {len(combined_artifacts['failed_grades'])} failed grades")
    
    return combined_artifacts

def download_and_extract_single_artifact(s3_path: str, output_dir: Path, artifact_type: str, force_download: bool = False, temp_version_dir: str = None) -> bool:
    """Download and extract a single artifact."""
    # Parse agent/task/filename from S3 path
    # e.g., s3://bucket/v1/artifacts/runs/stella/polarishub-tdcommons-herg/file.tar.gz
    path_parts = s3_path.replace('s3://', '').split('/')
    
    # If temp_version_dir is specified, use it to separate v1/v2 downloads
    if temp_version_dir:
        type_dir = output_dir / temp_version_dir / artifact_type
    else:
        type_dir = output_dir / artifact_type
    
    if len(path_parts) >= 6:  # organized structure
        agent = path_parts[4]
        task_safe = path_parts[5]
        filename = path_parts[6]
        
        # Create agent/task directory structure locally
        local_dir = type_dir / agent / task_safe
        local_dir.mkdir(parents=True, exist_ok=True)
        local_file = local_dir / filename
        
        # Check if extracted content already exists (skip re-download unless forced)
        if not force_download:
            if filename.endswith('.tar.gz'):
                extracted_name = filename.replace('.tar.gz', '')
                if (local_dir / extracted_name).exists():
                    return True  # Already extracted
            elif filename.endswith('.gz'):
                uncompressed_name = filename.replace('.gz', '')
                if (local_dir / uncompressed_name).exists():
                    return True  # Already decompressed
    else:  # flat structure fallback
        filename = path_parts[-1]
        local_file = type_dir / filename
        
        # Check if extracted content already exists
        if not force_download:
            if filename.endswith('.tar.gz'):
                extracted_name = filename.replace('.tar.gz', '')
                if (type_dir / extracted_name).exists():
                    return True
            elif filename.endswith('.gz'):
                uncompressed_name = filename.replace('.gz', '')
                if (type_dir / uncompressed_name).exists():
                    return True
    
    # Skip download if file already exists locally (unless forced)
    if not force_download and local_file.exists():
        return True
    
    # Download file
    exit_code, output = run_command([
        "aws", "s3", "cp", s3_path, str(local_file)
    ])
    
    if exit_code != 0:
        raise RuntimeError(f"Failed to download {s3_path}: {output}")
    
    # Extract if compressed
    if filename.endswith('.tar.gz'):
        with tarfile.open(local_file, 'r:gz') as tar:
            tar.extractall(local_file.parent, filter='data')
        local_file.unlink()  # Remove compressed file after extraction

    
    elif filename.endswith('.gz') and not filename.endswith('.tar.gz'):
        uncompressed_file = local_file.with_suffix('')
        with gzip.open(local_file, 'rb') as f_in:
            with open(uncompressed_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        local_file.unlink()  # Remove compressed file after extraction

    
    return True

def download_and_extract_artifacts(artifacts: Dict[str, List[str]], output_dir: Path, max_workers: int = 10, force_download: bool = False, v1v2_compare: bool = False):
    """Download and extract all artifacts in parallel.
    
    For failed_runs: Downloads v1 and v2 separately, then combines them properly
    to ensure v1 manual tasks are replaced by v2 manual tasks.
    """
    output_dir.mkdir(exist_ok=True)
    
    # Separate v1 and v2 artifacts
    v1_artifacts = {'runs': [], 'grades': [], 'failed_runs': [], 'failed_grades': []}
    v2_artifacts = {'runs': [], 'grades': [], 'failed_runs': [], 'failed_grades': []}
    
    for artifact_type, s3_paths in artifacts.items():
        for s3_path in s3_paths:
            if isinstance(s3_path, tuple):
                # Handle tuples from flat structure marking
                s3_path = s3_path[0]
            
            if '/v1/artifacts/' in s3_path:
                v1_artifacts[artifact_type].append(s3_path)
            elif '/v2/artifacts/' in s3_path:
                v2_artifacts[artifact_type].append(s3_path)
    
    # Download v1 artifacts
    log("Downloading v1 artifacts...")
    download_tasks = []
    for artifact_type, s3_paths in v1_artifacts.items():
        for s3_path in s3_paths:
            download_tasks.append((s3_path, output_dir, artifact_type, 'v1_temp' if artifact_type == 'failed_runs' else None))
    
    if download_tasks:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(download_and_extract_single_artifact, s3_path, output_dir, artifact_type, force_download, temp_dir): s3_path
                for s3_path, output_dir, artifact_type, temp_dir in download_tasks
            }
            for future in as_completed(futures):
                future.result()
    
    # Download v2 artifacts
    log("Downloading v2 artifacts...")
    download_tasks = []
    for artifact_type, s3_paths in v2_artifacts.items():
        for s3_path in s3_paths:
            download_tasks.append((s3_path, output_dir, artifact_type, 'v2_temp' if artifact_type == 'failed_runs' else None))
    
    if download_tasks:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(download_and_extract_single_artifact, s3_path, output_dir, artifact_type, force_download, temp_dir): s3_path
                for s3_path, output_dir, artifact_type, temp_dir in download_tasks
            }
            for future in as_completed(futures):
                future.result()
    
    # Now handle failed_runs: Remove manual tasks from v1, copy manual tasks from v2
    # Skip this special handling in v1v2_compare mode
    if v1v2_compare:
        # For v1v2 compare, just move the temp directories to final location
        for artifact_type in ['runs', 'grades', 'failed_runs', 'failed_grades']:
            temp_dir = output_dir / 'v1_temp' / artifact_type
            final_dir = output_dir / artifact_type
            if temp_dir.exists():
                if final_dir.exists():
                    shutil.rmtree(final_dir)
                shutil.move(str(temp_dir), str(final_dir))
        # Clean up temp directory
        temp_path = output_dir / 'v1_temp'
        if temp_path.exists():
            shutil.rmtree(temp_path)
        return
    
    v1_failed_runs_dir = output_dir / 'v1_temp' / 'failed_runs'
    v2_failed_runs_dir = output_dir / 'v2_temp' / 'failed_runs'
    final_failed_runs_dir = output_dir / 'failed_runs'
    
    if v1_failed_runs_dir.exists() or v2_failed_runs_dir.exists():
        log("Combining v1 and v2 failed runs (excluding v1 manual tasks)...")
        
        # First, move v1 failed runs to final location (if they exist)
        if v1_failed_runs_dir.exists():
            # Delete any existing failed_runs directory
            if final_failed_runs_dir.exists():
                shutil.rmtree(final_failed_runs_dir)
            
            # Move v1 to final location
            shutil.move(str(v1_failed_runs_dir), str(final_failed_runs_dir))
            
            # Now delete all manual tasks from v1
            log("Removing v1 manual task failures...")
            removed_count = 0
            for item in final_failed_runs_dir.rglob('*'):
                if item.is_dir() and '/manual/' in str(item):
                    # This is a manual task directory - delete it
                    shutil.rmtree(item)
                    removed_count += 1
                    log(f"  Removed: {item.name}")
            log(f"Removed {removed_count} v1 manual task failed runs")
        
        # Now copy v2 manual tasks to final location
        if v2_failed_runs_dir.exists():
            log("Copying v2 manual task failures...")
            copied_count = 0
            for item in v2_failed_runs_dir.rglob('*'):
                if item.is_dir() and '/manual/' in str(item):
                    # This is a manual task directory - copy it
                    relative_path = item.relative_to(v2_failed_runs_dir)
                    dest_path = final_failed_runs_dir / relative_path
                    
                    # Ensure parent directory exists
                    dest_path.parent.mkdir(parents=True, exist_ok=True)
                    
                    # Copy the directory
                    if dest_path.exists():
                        shutil.rmtree(dest_path)
                    shutil.copytree(str(item), str(dest_path))
                    copied_count += 1
                    log(f"  Copied: {item.name}")
            log(f"Copied {copied_count} v2 manual task failed runs")
        
        # Clean up temp directories
        if v1_failed_runs_dir.parent.exists():
            shutil.rmtree(v1_failed_runs_dir.parent)
        if v2_failed_runs_dir.parent.exists():
            shutil.rmtree(v2_failed_runs_dir.parent)
    
    log("✅ Download and merge complete!")

def count_total_run_attempts(runs_dir: Path, failed_runs_dir: Path, task_filter: List[str] = None) -> Dict[str, Dict[str, int]]:
    """Count total run attempts per agent-task from both successful and failed runs.
    
    Returns:
        Dict with run counts AND set of excluded run IDs
    """
    log("Counting total run attempts from runs and failed_runs directories...")
    
    run_counts = {}  # {agent_task_key: {'total': int, 'successful': int, 'failed': int}}
    excluded_quota_failures = 0
    excluded_quick_failures = 0
    excluded_max_steps = 0
    excluded_grading_server_failures = 0
    excluded_run_ids = set()  # Track which agent-task-run combinations were excluded
    
    # Process successful runs
    if runs_dir.exists():
        for metadata_file in runs_dir.rglob('metadata.json'):
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)
            
            if 'runs' not in metadata:
                continue
            
            for run_id, run_data in metadata['runs'].items():
                if 'task_id' not in run_data or 'agent_id' not in run_data:
                    continue
                
                task_id = run_data['task_id']
                agent_id = run_data['agent_id']
                
                # Apply task filter if specified
                if task_filter:
                    task_matches = any(pattern in task_id for pattern in task_filter)
                    if not task_matches:
                        continue
                
                # Check for quick failure even in successful runs (might have completed but with system error)
                # Look for run.log in the actual run directory
                run_dir = metadata_file.parent
                run_log_path = None
                
                # Try to find run.log in the current directory or subdirectories
                for potential_log in run_dir.rglob('run.log'):
                    run_log_path = potential_log
                    break
                
                if run_log_path:
                    if check_for_quick_failure(run_log_path, agent_id):
                        excluded_quick_failures += 1
                        # Track this excluded run by agent-task-run combination
                        run_dir_name = str(metadata_file.parent.name)
                        exclusion_key = f"{agent_id}::{task_id}::{run_dir_name}"
                        excluded_run_ids.add(exclusion_key)
                        continue  # Skip this run due to quick failure
                    
                    if check_for_max_steps_error(run_log_path, agent_id):
                        excluded_max_steps += 1
                        # Track this excluded run by agent-task-run combination
                        run_dir_name = str(metadata_file.parent.name)
                        exclusion_key = f"{agent_id}::{task_id}::{run_dir_name}"
                        excluded_run_ids.add(exclusion_key)
                        continue  # Skip this run due to max steps error
                    
                    if check_for_grading_server_failure(run_log_path):
                        excluded_grading_server_failures += 1
                        # Track this excluded run by agent-task-run combination
                        run_dir_name = str(metadata_file.parent.name)
                        exclusion_key = f"{agent_id}::{task_id}::{run_dir_name}"
                        excluded_run_ids.add(exclusion_key)
                        continue  # Skip this run due to grading server failure
                
                agent_task_key = f"{agent_id}::{task_id}"
                if agent_task_key not in run_counts:
                    run_counts[agent_task_key] = {'total': 0, 'successful': 0, 'failed': 0, 'agent': agent_id, 'task_id': task_id}
                
                run_counts[agent_task_key]['total'] += 1
                run_counts[agent_task_key]['successful'] += 1
                    

    
    # Process failed runs
    if failed_runs_dir.exists():
        for metadata_file in failed_runs_dir.rglob('metadata.json'):
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)
            
            if 'runs' not in metadata:
                continue
            
            for run_id, run_data in metadata['runs'].items():
                if 'task_id' not in run_data or 'agent_id' not in run_data:
                    continue
                
                task_id = run_data['task_id']
                agent_id = run_data['agent_id']
                
                # Apply task filter if specified
                if task_filter:
                    task_matches = any(pattern in task_id for pattern in task_filter)
                    if not task_matches:
                        continue
                
                # Check for API quota/credit errors in run.log
                # Look for run.log in the actual run directory
                run_dir = metadata_file.parent
                run_log_path = None
                
                # Try to find run.log in the current directory or subdirectories
                for potential_log in run_dir.rglob('run.log'):
                    run_log_path = potential_log
                    break
                
                if run_log_path:
                    if check_for_api_quota_errors(run_log_path):
                        excluded_quota_failures += 1
                        run_dir_name = str(metadata_file.parent.name)
                        exclusion_key = f"{agent_id}::{task_id}::{run_dir_name}"
                        excluded_run_ids.add(exclusion_key)
                        continue  # Skip this failed run due to quota/credit issues
                    
                    # Check for quick failure (< 300 seconds)
                    if check_for_quick_failure(run_log_path, agent_id):
                        excluded_quick_failures += 1
                        run_dir_name = str(metadata_file.parent.name)
                        exclusion_key = f"{agent_id}::{task_id}::{run_dir_name}"
                        excluded_run_ids.add(exclusion_key)
                        continue  # Skip this failed run due to quick failure
                    
                    # Check for max steps error (STELLA only)
                    if check_for_max_steps_error(run_log_path, agent_id):
                        excluded_max_steps += 1
                        run_dir_name = str(metadata_file.parent.name)
                        exclusion_key = f"{agent_id}::{task_id}::{run_dir_name}"
                        excluded_run_ids.add(exclusion_key)
                        continue  # Skip this failed run due to max steps error
                    
                    # Check for grading server failure
                    if check_for_grading_server_failure(run_log_path):
                        excluded_grading_server_failures += 1
                        run_dir_name = str(metadata_file.parent.name)
                        exclusion_key = f"{agent_id}::{task_id}::{run_dir_name}"
                        excluded_run_ids.add(exclusion_key)
                        continue  # Skip this failed run due to grading server failure
                else:
                    continue
                
                agent_task_key = f"{agent_id}::{task_id}"
                if agent_task_key not in run_counts:
                    run_counts[agent_task_key] = {'total': 0, 'successful': 0, 'failed': 0, 'agent': agent_id, 'task_id': task_id}
                
                run_counts[agent_task_key]['total'] += 1
                run_counts[agent_task_key]['failed'] += 1
                

    
    log(f"Found run attempts for {len(run_counts)} agent-task combinations")
    if excluded_quota_failures > 0:
        log(f"🚫 Excluded {excluded_quota_failures} failed runs due to API quota/credit issues")
    if excluded_quick_failures > 0:
        log(f"🚫 Excluded {excluded_quick_failures} runs due to quick failure (< 300s runtime, except dummy agent)")
    if excluded_max_steps > 0:
        log(f"🚫 Excluded {excluded_max_steps} STELLA runs due to max steps error")
    if excluded_grading_server_failures > 0:
        log(f"🚫 Excluded {excluded_grading_server_failures} runs due to grading server failure")
    
    return run_counts, excluded_run_ids

def calculate_ranks_for_results(all_results: List[Dict]) -> List[Dict]:
    """Calculate ranks for each task based on average leaderboard percentile.
    
    For each task:
    1. Calculate average leaderboard percentile per agent
    2. Rank agents (1=best/highest percentile, n=worst/lowest)
    3. Handle ties with fractional ranks
    """
    log("Calculating ranks for each task...")
    
    # Convert to DataFrame for easier manipulation
    df = pd.DataFrame(all_results)
    
    # Calculate average percentile per agent-task
    task_avg = df.groupby(['agent', 'task_id'])['leaderboard_percentile'].mean().reset_index()
    
    # Calculate ranks for each task
    ranks_data = []
    for task_id in task_avg['task_id'].unique():
        task_data = task_avg[task_avg['task_id'] == task_id].copy()
        
        # Rank agents for this task (higher percentile = better = lower rank number)
        # Use 'average' method for ties (fractional ranks)
        # Negate percentiles so that higher percentile gets lower (better) rank
        task_data['rank'] = rankdata(-task_data['leaderboard_percentile'], method='average')
        
        ranks_data.append(task_data)
    
    # Combine all rank data
    ranks_df = pd.concat(ranks_data, ignore_index=True)
    
    # Merge ranks back to all results
    df = df.merge(ranks_df[['agent', 'task_id', 'rank']], on=['agent', 'task_id'], how='left')
    
    # Convert back to list of dicts
    results_with_ranks = df.to_dict('records')
    
    # Log summary statistics
    rank_summary = df.groupby('agent')['rank'].agg(['mean', 'std']).round(2)
    log("Average ranks by agent:")
    for agent, row in rank_summary.iterrows():
        log(f"  {agent}: {row['mean']:.2f} ± {row['std']:.2f}")
    
    return results_with_ranks

def analyze_grading_results(grades_dir: Path, runs_dir: Path, failed_runs_dir: Path, task_filter: List[str] = None) -> Dict[str, Any]:
    """Extract ALL individual results and compute completion rates."""
    log("Extracting individual grading results...")
    
    if task_filter:
        log(f"Applying task filter: {task_filter}")
    
    # First, count total run attempts from both successful and failed runs
    run_counts, excluded_run_ids = count_total_run_attempts(runs_dir, failed_runs_dir, task_filter)
    
    # Store ALL individual results (replicates) - NO AGGREGATION
    all_results = []
    # Track successful completions with valid scores by agent-task combination
    successful_completions = {}
    
    # Find all individual report JSON files (extracted from individual_reports directories)
    report_files = []
    original_label_projection_files = []
    regraded_label_projection_files = []
    
    # Include original grading results
    for root, dirs, files in os.walk(grades_dir):
        for file in files:
            # Look for individual task report JSON files (extracted from individual_reports.tar.gz)
            if file.endswith('.json') and not file.endswith('_grading_report.json'):
                file_path = Path(root) / file
                
                # Check if this is a label-projection result
                with open(file_path, 'r') as f:
                    report_data = json.load(f)
                
                if report_data.get('task_id') == 'manual/open-problems-label-projection':
                    original_label_projection_files.append(file_path)
                else:
                    report_files.append(file_path)
    
    # Include regraded label-projection results (these replace original ones)
    regraded_dir = grades_dir.parent / 'regraded_label_projection'
    if regraded_dir.exists():
        log(f"Including regraded label-projection results from {regraded_dir}")
        # Group regraded results by run identifier
        regraded_by_run = {}
        for root, dirs, files in os.walk(regraded_dir):
            for file in files:
                # Look for individual task report JSON files from regrading
                if file.endswith('.json') and not file.endswith('_grading_report.json') and 'individual_reports' in root:
                    file_path = Path(root) / file
                    # Extract run identifier (e.g., 2025-08-23T14-42-13-GMT_run-group_aide)
                    run_id = file_path.parts[-4]  # Get the run directory name
                    # Extract timestamp of this regrading (e.g., 2025-08-24T15-35-39-GMT)
                    regrade_timestamp = file_path.parts[-2].replace('_individual_reports', '')
                    
                    if run_id not in regraded_by_run:
                        regraded_by_run[run_id] = []
                    regraded_by_run[run_id].append((regrade_timestamp, file_path))
        
        # Use only the most recent regrading for each run
        for run_id, regrade_list in regraded_by_run.items():
            # Sort by timestamp and take the most recent
            regrade_list.sort(key=lambda x: x[0], reverse=True)
            most_recent_path = regrade_list[0][1]
            regraded_label_projection_files.append(most_recent_path)
            if len(regrade_list) > 1:
                log(f"  {run_id}: Using most recent regrading from {regrade_list[0][0]} (ignoring {len(regrade_list)-1} older regradings)")
    
    # Use regraded results if available, otherwise use original
    if regraded_label_projection_files:
        log(f"🔄 Using {len(regraded_label_projection_files)} regraded label-projection results instead of {len(original_label_projection_files)} original ones")
        report_files.extend(regraded_label_projection_files)
    else:
        log(f"📊 Using {len(original_label_projection_files)} original label-projection results")
        report_files.extend(original_label_projection_files)
    
    log(f"Found {len(report_files)} individual result files to process")
    
    for report_file in report_files:
        with open(report_file, 'r') as f:
            report = json.load(f)
        
        # Extract task_id early for exclusion check
        if 'task_id' not in report:
            raise ValueError(f"Missing required field 'task_id' in {report_file}")
        task_id = report['task_id']
        
        # Extract agent from file path
        # e.g., grades/stella/task_name/... -> stella
        agent = None
        path_parts = Path(report_file).parts
        if 'grades' in path_parts:
            grades_idx = path_parts.index('grades')
            if grades_idx + 1 < len(path_parts):
                agent = path_parts[grades_idx + 1]
        
        # Check if this grading result is from an excluded run using the submission_path field
        if 'submission_path' in report and report['submission_path']:
            submission_path = report['submission_path']
            # Extract the run ID from the submission path
            # Format: /home/runner/biomlbench/runs/2025-08-20T04-23-40-GMT_run-group_stella/...
            run_id_match = None
            for part in submission_path.split('/'):
                if '_run-group_' in part:
                    run_id_match = part
                    break
            
            if run_id_match and agent and task_id:
                # Check if this specific agent-task-run combination was excluded
                exclusion_key = f"{agent}::{task_id}::{run_id_match}"
                if exclusion_key in excluded_run_ids:
                    log(f"⚠️  Skipping grading result from excluded run: {exclusion_key}")
                    continue
        
        # REQUIRED fields - no defaults
        if 'score' not in report:
            raise ValueError(f"Missing required field 'score' in {report_file}")
        if 'leaderboard_percentile' not in report:
            raise ValueError(f"Missing required field 'leaderboard_percentile' in {report_file}")
        if 'gold_medal' not in report:
            raise ValueError(f"Missing required field 'gold_medal' in {report_file}")
        if 'silver_medal' not in report:
            raise ValueError(f"Missing required field 'silver_medal' in {report_file}")
        if 'bronze_medal' not in report:
            raise ValueError(f"Missing required field 'bronze_medal' in {report_file}")
        if 'above_median' not in report:
            raise ValueError(f"Missing required field 'above_median' in {report_file}")
        
        # task_id and agent already extracted above
        score = report['score']
        
        # Check if this is from regraded results: .../regraded_label_projection/2025-XX-XXTXX-XX-XX-GMT_run-group_AGENT/...
        if agent is None:
            path_parts = report_file.parts
            for i, part in enumerate(path_parts):
                if part == 'regraded_label_projection' and i + 1 < len(path_parts):
                    run_group_part = path_parts[i + 1]
                    if '_run-group_' in run_group_part:
                        agent = run_group_part.split('_run-group_')[-1]
                        break
        
        if agent is None:
            raise ValueError(f"Could not extract agent from file path in {report_file}")
        
        # Apply task filter if specified
        if task_filter:
            task_matches = any(pattern in task_id for pattern in task_filter)
            if not task_matches:
                continue  # Skip this task if it doesn't match filter
        
        # Track successful completions with valid scores
        agent_task_key = f"{agent}::{task_id}"
        if agent_task_key not in successful_completions:
            successful_completions[agent_task_key] = set()  # Use set to avoid double-counting
        
        if score is None:
            log(f"⚠️  Score is null in {report_file} - treating as 0")
            score = 0.0
            # For null scores, set metrics appropriately
            percentile = 0.0  # Null score gets 0 percentile
            gold_medal = False
            silver_medal = False
            bronze_medal = False
            above_median = False
        else:
            # Count as successful completion (only if score is not null)
            # Extract run ID from submission_path to avoid double-counting multiple gradings of the same run
            if 'submission_path' in report and report['submission_path']:
                # Extract run ID from submission path (e.g., 2025-08-20T04-23-40-GMT_run-group_stella)
                submission_path = report['submission_path']
                run_id_match = None
                for part in submission_path.split('/'):
                    if '_run-group_' in part:
                        run_id_match = part
                        break
                if run_id_match:
                    successful_completions[agent_task_key].add(run_id_match)
                else:
                    log(f"⚠️  Could not extract run ID from submission_path: {submission_path}")
            else:
                log(f"⚠️  No submission_path in grading report: {report_file}")
            
            # Convert types explicitly for non-null scores
            percentile = float(report['leaderboard_percentile'])
            gold_medal = bool(report['gold_medal'])
            silver_medal = bool(report['silver_medal'])
            bronze_medal = bool(report['bronze_medal'])
            above_median = bool(report['above_median'])
        
        # Store this individual result - each JSON file is one replicate
        result = {
            'agent': agent,
            'task_id': task_id,
            'score': float(score),
            'leaderboard_percentile': percentile,
            'above_median': above_median,
            'gold_medal': gold_medal,
            'silver_medal': silver_medal,
            'bronze_medal': bronze_medal,
            'any_medal': gold_medal or silver_medal or bronze_medal,
            'report_file': str(report_file)
        }
        all_results.append(result)
            
    log(f"Extracted {len(all_results)} individual results")
    
    # Calculate completion rates using total run attempts and successful completions
    completion_rates = []
    for agent_task_key, counts in run_counts.items():
        successful_runs = successful_completions.get(agent_task_key, set())
        successful_count = len(successful_runs)  # Count unique successful runs
        total_attempts = counts['total']
        completion_rate = successful_count / total_attempts if total_attempts > 0 else 0.0
        
        completion_rates.append({
            'agent': counts['agent'],
            'task_id': counts['task_id'],
            'total_attempts': total_attempts,
            'successful_completions': successful_count,
            'completion_rate': completion_rate
        })
    
    log(f"Calculated completion rates for {len(completion_rates)} agent-task combinations")
    
    # Add zero-score entries for all failed attempts without grading reports
    log("Adding zero-score entries for failed attempts without grading reports...")
    added_zeros = 0
    for agent_task_key, counts in run_counts.items():
        successful_runs = successful_completions.get(agent_task_key, set())
        successful_count = len(successful_runs)  # Count unique successful runs
        failed_count = counts['total'] - successful_count
        
        if failed_count > 0:
            # Add zero-score entries for each failed attempt
            for _ in range(failed_count):
                all_results.append({
                    'agent': counts['agent'],
                    'task_id': counts['task_id'],
                    'score': 0.0,
                    'leaderboard_percentile': 0.0,
                    'above_median': False,
                    'gold_medal': False,
                    'silver_medal': False,
                    'bronze_medal': False,
                    'any_medal': False,
                    'report_file': 'failed_run_no_report'
                })
                added_zeros += 1
    
    log(f"Added {added_zeros} zero-score entries for failed runs")
    
    # Calculate ranks for all results
    all_results = calculate_ranks_for_results(all_results)
    
    # Show completion summary
    if completion_rates:
        total_attempts = sum(cr['total_attempts'] for cr in completion_rates)
        total_successes = sum(cr['successful_completions'] for cr in completion_rates)
        overall_completion_rate = total_successes / total_attempts if total_attempts > 0 else 0.0
        log(f"Overall completion rate: {total_successes}/{total_attempts} ({overall_completion_rate:.1%})")
    
    return {
        'all_results': all_results,
        'completion_rates': completion_rates
    }

def analyze_run_metadata(runs_dir: Path) -> Dict[str, Any]:
    """Analyze run metadata for execution times and resource usage."""
    log("Analyzing run metadata...")
    
    metadata_files = []
    for root, dirs, files in os.walk(runs_dir):
        for file in files:
            if file == 'metadata.json':
                metadata_files.append(Path(root) / file)
    
    analysis = {
        'execution_times': {},
        'task_durations': {},
        'agent_performance': {}
    }
    
    for metadata_file in metadata_files:
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)
        
        # Require essential fields
        if 'agent_id' not in metadata:
            raise ValueError(f"Missing required field 'agent_id' in {metadata_file}")
        if 'runs' not in metadata:
            raise ValueError(f"Missing required field 'runs' in {metadata_file}")
        
        agent_id = metadata['agent_id']
        runs = metadata['runs']
        
        # Initialize agent performance tracking
        if agent_id not in analysis['agent_performance']:
            analysis['agent_performance'][agent_id] = {
                'total_time': 0.0,
                'successful_runs': 0,
                'failed_runs': 0
            }
        
        for run_id, run_data in runs.items():
            if 'task_id' not in run_data:
                raise ValueError(f"Missing required field 'task_id' in run {run_id} of {metadata_file}")
            if 'success' not in run_data:
                raise ValueError(f"Missing required field 'success' in run {run_id} of {metadata_file}")
            
            task_id = run_data['task_id']
            success = run_data['success']
            
            if success:
                analysis['agent_performance'][agent_id]['successful_runs'] += 1
            else:
                analysis['agent_performance'][agent_id]['failed_runs'] += 1
            

    return analysis

# Removed generate_summary_report - only raw CSV output needed

def generate_csv_reports(all_results: List[Dict], completion_rates: List[Dict], output_dir: Path):
    """Generate CSV from raw results and completion rates."""
    log("Generating CSV reports...")
    
    # Replicate-level results CSV - this is the PRIMARY output
    if len(all_results) > 0:
        replicate_df = pd.DataFrame(all_results)
        
        # Remove internal fields before saving
        if 'report_file' in replicate_df.columns:
            replicate_df = replicate_df.drop('report_file', axis=1)
        
        replicate_df.to_csv(output_dir / 'all_replicates.csv', index=False)
        log(f"Saved {len(all_results)} individual results to all_replicates.csv")
    else:
        log("⚠️  No successful results to save to all_replicates.csv")
    
    # Completion rates CSV - NEW output
    if len(completion_rates) > 0:
        completion_df = pd.DataFrame(completion_rates)
        completion_df.to_csv(output_dir / 'completion_rates.csv', index=False)
        log(f"Saved completion rates for {len(completion_rates)} agent-task combinations to completion_rates.csv")
    else:
        log("⚠️  No completion rate data to save")

def main():
    parser = argparse.ArgumentParser(
        description="Download and analyze BioML-bench results from S3"
    )
    parser.add_argument(
        "--bucket",
        default="biomlbench",
        help="S3 bucket name (default: biomlbench)"
    )
    parser.add_argument(
        "--prefix", 
        default="v1/artifacts",
        help="S3 prefix (default: v1/artifacts)"
    )
    parser.add_argument(
        "--output-dir",
        default="analysis_results",
        help="Local directory to store downloaded artifacts and analysis (default: analysis_results)"
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Skip download phase and analyze existing local data"
    )
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Force re-download even if files already exist locally"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true", 
        help="Show what would be downloaded without actually downloading"
    )
    parser.add_argument(
        "--parallel-downloads",
        type=int,
        default=10,
        help="Number of parallel download workers (default: 10)"
    )
    parser.add_argument(
        "--agents",
        nargs='+',
        help="Filter to specific agents (e.g., --agents stella mlagentbench)"
    )
    parser.add_argument(
        "--tasks",
        nargs='+', 
        help="Filter to specific task patterns (e.g., --tasks polarishub kaggle)"
    )
    parser.add_argument(
        "--fix-label-projection",
        action="store_true",
        help="Fix submission paths and regrade label-projection tasks"
    )
    parser.add_argument(
        "--v1v2-compare",
        action="store_true",
        help="Use v1 or v2 artifacts only (no merging) - use with --prefix v1/artifacts or v2/artifacts"
    )
    
    args = parser.parse_args()
    output_dir = Path(args.output_dir)
    
    if not args.skip_download:
        # Discover artifacts
        if args.v1v2_compare:
            # Use single version discovery for v1v2 comparison mode
            artifacts = discover_s3_artifacts_single_version(args.bucket, args.prefix, args.agents, args.tasks, exclude_tasks=[])
        else:
            artifacts = discover_s3_artifacts(args.bucket, args.prefix, args.agents, args.tasks)
        
        if args.dry_run:
            log("DRY RUN - Would download:")
            for artifact_type, paths in artifacts.items():
                log(f"  {artifact_type}: {len(paths)} files")
                for path in paths[:3]:  # Show first 3 examples
                    log(f"    {path}")
                if len(paths) > 3:
                    log(f"    ... and {len(paths) - 3} more")
            return 0
        
        # Download and extract
        download_and_extract_artifacts(artifacts, output_dir, args.parallel_downloads, args.force_download, args.v1v2_compare)
    
    # Fix label projection submissions if requested
    if args.fix_label_projection:
        runs_dir = output_dir / 'runs'
        if runs_dir.exists():
            fix_submission_paths_and_regrade(runs_dir, output_dir)
        else:
            log("⚠️  No runs directory found for fixing label-projection submissions")
    
    # Analyze results
    grades_dir = output_dir / 'grades'
    runs_dir = output_dir / 'runs'
    failed_runs_dir = output_dir / 'failed_runs'
    
    if grades_dir.exists():
        results_data = analyze_grading_results(grades_dir, runs_dir, failed_runs_dir, args.tasks)
        all_results = results_data['all_results']
        completion_rates = results_data['completion_rates']
        
        # Generate CSV with ALL individual results AND completion rates
        generate_csv_reports(all_results, completion_rates, output_dir)
        
        log("Analysis complete!")
        log(f"📈 Individual results CSV: {output_dir / 'all_replicates.csv'}")
        log(f"📊 Completion rates CSV: {output_dir / 'completion_rates.csv'}")
        
        # Quick terminal summary
        print("\n" + "="*60)
        print("RESULTS SUMMARY")
        print("="*60)
        
        if len(all_results) == 0:
            raise ValueError("No results found!")
        
        total_replicates = len(all_results)
        unique_agents = set(r['agent'] for r in all_results)
        unique_tasks = set(r['task_id'] for r in all_results)
        medal_winners = sum(1 for r in all_results if r['any_medal'])
        
        print(f"Total replicates: {total_replicates}")
        print(f"Unique agents: {len(unique_agents)} ({', '.join(sorted(unique_agents))})")
        print(f"Unique tasks: {len(unique_tasks)}")
        print(f"Medal-winning replicates: {medal_winners}")
        
        # Show replicate counts per agent-task
        from collections import Counter
        agent_task_counts = Counter((r['agent'], r['task_id']) for r in all_results)
        print(f"\nReplicate counts per agent-task:")
        for (agent, task), count in sorted(agent_task_counts.items()):
            print(f"  {agent} x {task}: {count} replicates")
        
        print("="*60)
        
    else:
        raise ValueError("No grading results directory found!")
    
    return 0

if __name__ == "__main__":
    exit(main()) 