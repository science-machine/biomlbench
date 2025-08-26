#!/usr/bin/env python3
"""
Fix All Submission Paths Script

Fixes submission paths in ALL submission.jsonl files by updating them to point to 
the correct local paths instead of the original runner environment paths.

This fixes completion rate discrepancies caused by grading failures due to path mismatches.
"""

import json
import subprocess
from pathlib import Path
import argparse
from typing import List, Dict

def log(message: str):
    """Simple logging with timestamp."""
    import time
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def run_command(cmd: List[str], description: str = "") -> tuple:
    """Run a command and return (exit_code, output)."""
    if description:
        log(f"Running: {description}")
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode, result.stdout + result.stderr

def fix_all_submission_paths_and_regrade(runs_dir: Path, output_dir: Path):
    """Fix submission paths in ALL submission.jsonl files and regrade them."""
    log("🔧 Fixing submission paths for ALL tasks...")
    
    # Find all submission.jsonl files
    all_submissions = []
    for submission_file in runs_dir.rglob('submission.jsonl'):
        with open(submission_file, 'r') as f:
            submission_data = json.load(f)
        all_submissions.append((submission_file, submission_data))
    
    log(f"Found {len(all_submissions)} submission.jsonl files to fix")
    
    if not all_submissions:
        return
    
    fixed_count = 0
    regrade_paths = []
    
    for submission_file, submission_data in all_submissions:
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
        
        # Hardcoded exclusion for problematic tasks
        task_id = submission_data.get('task_id', '')
        if task_id == 'polarishub/polaris-adme-fang-r-1':
            continue
        
        # Check if paths need fixing (contain '/home/runner/biomlbench')
        old_submission_path = submission_data.get('submission_path', '')
        if '/home/runner/biomlbench' not in old_submission_path:
            # Paths are already correct
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
        log(f"🎯 Regrading {len(regrade_paths)} fixed submissions...")
        
        for submission_file in regrade_paths:
            log(f"Regrading {submission_file}")
            
            # Create output directory for this grading
            # Use the task name and run timestamp for organization
            with open(submission_file, 'r') as f:
                sub_data = json.load(f)
            
            task_id = sub_data['task_id'].replace('/', '-')
            run_timestamp = submission_file.parent.name
            
            grade_output_dir = output_dir / 'regraded_all_tasks' / task_id / run_timestamp
            grade_output_dir.mkdir(parents=True, exist_ok=True)
            
            # Run biomlbench grade command
            cmd = [
                'uv', 'run', 'biomlbench', 'grade',
                '--submission', str(submission_file),
                '--output-dir', str(grade_output_dir)
            ]
            
            exit_code, output = run_command(cmd, f"Regrade {submission_file.name}")
            
            if exit_code == 0:
                log(f"✅ Successfully regraded {submission_file.name}")
            else:
                log(f"❌ Failed to regrade {submission_file.name}: {output}")
        
        log("🎉 Regrading complete!")
        
        # Show breakdown by task
        task_breakdown = {}
        for submission_file in regrade_paths:
            with open(submission_file, 'r') as f:
                sub_data = json.load(f)
            task_id = sub_data['task_id']
            if task_id not in task_breakdown:
                task_breakdown[task_id] = 0
            task_breakdown[task_id] += 1
        
        log("📈 Fixed submissions by task:")
        for task, count in sorted(task_breakdown.items()):
            log(f"  {task}: {count} submissions")

def main():
    parser = argparse.ArgumentParser(
        description="Fix submission paths for ALL tasks and regrade them"
    )
    parser.add_argument(
        "--analysis-dir",
        default="analysis_results_v1v2_combined",
        help="Analysis directory containing runs/ (default: analysis_results_v1v2_combined)"
    )
    
    args = parser.parse_args()
    
    analysis_dir = Path(args.analysis_dir)
    if not analysis_dir.exists():
        log(f"❌ Analysis directory does not exist: {analysis_dir}")
        return 1
    
    runs_dir = analysis_dir / 'runs'
    if not runs_dir.exists():
        log(f"❌ Runs directory does not exist: {runs_dir}")
        return 1
    
    fix_all_submission_paths_and_regrade(runs_dir, analysis_dir)
    
    return 0

if __name__ == "__main__":
    exit(main()) 