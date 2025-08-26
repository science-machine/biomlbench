#!/usr/bin/env python3
"""
Find Failed Grading Script

Identifies cases where:
1. Agent completed execution successfully (run exists in runs/ directory)
2. But grading failed (no successful grading report in grades/ directory)

This helps explain completion rates < 1.0 where failed runs aren't in failed_runs/ folder.
"""

import json
import shutil
from pathlib import Path
import argparse
from typing import List, Dict, Set

def log(message: str):
    """Simple logging with timestamp."""
    import time
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def extract_run_id_from_path(path_str: str) -> str:
    """Extract run ID from a file path."""
    # Look for pattern like 2025-08-23T04-48-11-GMT_run-group_biomni
    parts = path_str.split('/')
    for part in parts:
        if '_run-group_' in part:
            return part
    return None

def find_successful_gradings(grades_dir: Path) -> Set[str]:
    """Find all runs that have successful grading reports."""
    log(f"Scanning for successful grading reports in {grades_dir}")
    
    successful_run_ids = set()
    
    if not grades_dir.exists():
        log("⚠️  Grades directory does not exist")
        return successful_run_ids
    
    # Include original grading results
    for json_file in grades_dir.rglob('*.json'):
        if 'individual_reports' in str(json_file) and not json_file.name.endswith('_grading_report.json'):
            try:
                with open(json_file, 'r') as f:
                    report = json.load(f)
                
                # Extract run ID from submission_path if available
                if 'submission_path' in report and report['submission_path']:
                    run_id = extract_run_id_from_path(report['submission_path'])
                    if run_id:
                        successful_run_ids.add(run_id)
            except Exception as e:
                log(f"⚠️  Error reading grading report {json_file}: {e}")
    
    # Include regraded results
    regraded_dir = grades_dir.parent / 'regraded_label_projection'
    if regraded_dir.exists():
        log(f"Including regraded results from {regraded_dir}")
        for json_file in regraded_dir.rglob('*.json'):
            if 'individual_reports' in str(json_file) and not json_file.name.endswith('_grading_report.json'):
                try:
                    with open(json_file, 'r') as f:
                        report = json.load(f)
                    
                    # Extract run ID from submission_path if available
                    if 'submission_path' in report and report['submission_path']:
                        run_id = extract_run_id_from_path(report['submission_path'])
                        if run_id:
                            successful_run_ids.add(run_id)
                except Exception as e:
                    log(f"⚠️  Error reading regraded report {json_file}: {e}")
    
    log(f"Found {len(successful_run_ids)} runs with successful grading reports")
    return successful_run_ids

def find_failed_gradings(runs_dir: Path, grades_dir: Path, output_dir: Path):
    """Find runs where agent completed but grading failed."""
    log("🔍 Finding cases of successful runs with failed grading...")
    
    if not runs_dir.exists():
        log("⚠️  Runs directory does not exist")
        return
    
    # Get all successful grading run IDs
    successful_grading_run_ids = find_successful_gradings(grades_dir)
    
    failed_grading_cases = []
    total_runs_found = 0
    
    # Scan all runs
    for metadata_file in runs_dir.rglob('metadata.json'):
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)
        
        if 'runs' not in metadata:
            continue
        
        agent_id = metadata.get('agent_id', 'unknown')
        
        for run_id, run_data in metadata['runs'].items():
            total_runs_found += 1
            
            task_id = run_data.get('task_id', 'unknown')
            success = run_data.get('success', False)
            
            # Extract run directory name from metadata file path
            run_dir_name = metadata_file.parent.name
            
            # Check if this run has successful grading
            if run_dir_name not in successful_grading_run_ids:
                # This is a case of successful run but failed grading
                failed_grading_cases.append({
                    'agent_id': agent_id,
                    'task_id': task_id,
                    'run_id': run_id,
                    'run_dir_name': run_dir_name,
                    'success': success,
                    'metadata_file': metadata_file,
                    'run_directory': metadata_file.parent
                })
                
                log(f"❌ Failed grading: {agent_id} | {task_id} | {run_dir_name}")
    
    log(f"📊 Summary:")
    log(f"  Total runs found: {total_runs_found}")
    log(f"  Runs with successful grading: {len(successful_grading_run_ids)}")
    log(f"  Runs with failed grading: {len(failed_grading_cases)}")
    
    if failed_grading_cases:
        # Create output directory
        failed_grading_dir = output_dir / 'failed_grading_cases'
        failed_grading_dir.mkdir(parents=True, exist_ok=True)
        
        # Copy failed grading cases
        log(f"📁 Copying {len(failed_grading_cases)} failed grading cases to {failed_grading_dir}")
        
        copied_count = 0
        for case in failed_grading_cases:
            # Create subdirectory structure: agent/task/run_dir_name
            agent_task_dir = failed_grading_dir / case['agent_id'] / case['task_id'].replace('/', '-')
            agent_task_dir.mkdir(parents=True, exist_ok=True)
            
            # Copy the entire run directory
            src_dir = case['run_directory']
            dst_dir = agent_task_dir / case['run_dir_name']
            
            if not dst_dir.exists():
                shutil.copytree(src_dir, dst_dir)
                copied_count += 1
        
        log(f"✅ Copied {copied_count} failed grading cases")
        
        # Generate summary report
        summary_file = failed_grading_dir / 'summary.json'
        with open(summary_file, 'w') as f:
            summary_data = {
                'total_runs_found': total_runs_found,
                'successful_gradings': len(successful_grading_run_ids),
                'failed_gradings': len(failed_grading_cases),
                'failed_cases': [
                    {
                        'agent_id': case['agent_id'],
                        'task_id': case['task_id'],
                        'run_dir_name': case['run_dir_name'],
                        'success': case['success']
                    }
                    for case in failed_grading_cases
                ]
            }
            json.dump(summary_data, f, indent=2)
        
        log(f"📄 Summary report saved to {summary_file}")
        
        # Show breakdown by agent and task
        agent_task_breakdown = {}
        for case in failed_grading_cases:
            key = f"{case['agent_id']}::{case['task_id']}"
            if key not in agent_task_breakdown:
                agent_task_breakdown[key] = 0
            agent_task_breakdown[key] += 1
        
        log("📈 Breakdown by agent-task:")
        for key, count in sorted(agent_task_breakdown.items()):
            agent, task = key.split('::')
            log(f"  {agent} | {task}: {count} failed gradings")
    
    else:
        log("✅ No failed grading cases found!")

def main():
    parser = argparse.ArgumentParser(
        description="Find cases where agent completed execution but grading failed"
    )
    parser.add_argument(
        "--analysis-dir",
        default="analysis_results_v1v2_combined",
        help="Analysis directory containing runs/ and grades/ (default: analysis_results_v1v2_combined)"
    )
    parser.add_argument(
        "--output-dir",
        help="Output directory for failed grading cases (default: <analysis-dir>/failed_grading_analysis)"
    )
    
    args = parser.parse_args()
    
    analysis_dir = Path(args.analysis_dir)
    if not analysis_dir.exists():
        log(f"❌ Analysis directory does not exist: {analysis_dir}")
        return 1
    
    runs_dir = analysis_dir / 'runs'
    grades_dir = analysis_dir / 'grades'
    
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = analysis_dir / 'failed_grading_analysis'
    
    find_failed_gradings(runs_dir, grades_dir, output_dir)
    
    return 0

if __name__ == "__main__":
    exit(main()) 