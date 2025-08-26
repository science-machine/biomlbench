#!/usr/bin/env python3
"""
Simple script to find runs that exist but don't have grading reports.
"""

import json
from pathlib import Path
import sys

def find_missing_grades(runs_dir: Path, grades_dir: Path):
    """Find runs that exist but don't have grading reports."""
    
    print(f"Checking for runs without grades...")
    print(f"Runs directory: {runs_dir}")
    print(f"Grades directory: {grades_dir}")
    print()
    
    missing_grades = []
    
    # Find all runs
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
            
            # Convert task_id to safe directory name (replace / with -)
            task_safe = task_id.replace('/', '-')
            
            # Check if grading report exists
            # Look for individual reports in grades/agent/task/
            grade_pattern = grades_dir / agent_id / task_safe
            
            # Look for any JSON files in this directory (excluding _grading_report.json)
            grading_found = False
            if grade_pattern.exists():
                for json_file in grade_pattern.rglob('*.json'):
                    if not json_file.name.endswith('_grading_report.json'):
                        # Check if this JSON file corresponds to this run
                        with open(json_file, 'r') as f:
                            try:
                                grade_data = json.load(f)
                                if 'submission_path' in grade_data:
                                    submission_path = grade_data['submission_path']
                                    # Check if this grading matches this run
                                    run_dir_name = metadata_file.parent.name
                                    if run_dir_name in submission_path:
                                        grading_found = True
                                        break
                            except:
                                continue
            
            if not grading_found:
                missing_grades.append({
                    'agent': agent_id,
                    'task': task_id,
                    'run_dir': str(metadata_file.parent),
                    'run_id': run_id
                })
    
    return missing_grades

def main():
    if len(sys.argv) != 2:
        print("Usage: python find_missing_grades.py <analysis_directory>")
        sys.exit(1)
    
    analysis_dir = Path(sys.argv[1])
    runs_dir = analysis_dir / 'runs'
    grades_dir = analysis_dir / 'grades'
    
    if not runs_dir.exists():
        print(f"Error: Runs directory not found: {runs_dir}")
        sys.exit(1)
        
    if not grades_dir.exists():
        print(f"Error: Grades directory not found: {grades_dir}")
        sys.exit(1)
    
    missing = find_missing_grades(runs_dir, grades_dir)
    
    print(f"Found {len(missing)} runs without grading reports:")
    print()
    
    for item in missing:
        print(f"Agent: {item['agent']}")
        print(f"Task: {item['task']}")
        print(f"Run Directory: {item['run_dir']}")
        print(f"Run ID: {item['run_id']}")
        print("-" * 50)

if __name__ == "__main__":
    main() 