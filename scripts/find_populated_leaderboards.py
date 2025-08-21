#!/usr/bin/env python3
"""
Script to find all PolarisHub tasks with leaderboards containing at least 5 lines.

A leaderboard with at least 5 lines means it has:
- 1 header line (teamName,score,submissionDate)
- At least 4 data rows (team submissions)

This indicates the task has sufficient activity/participation.
"""

import os
import csv
from pathlib import Path

def count_non_empty_lines(file_path):
    """Count non-empty lines in a CSV file."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = [line.strip() for line in f.readlines()]
            # Filter out empty lines
            non_empty_lines = [line for line in lines if line]
            return len(non_empty_lines)
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return 0

def find_populated_leaderboards():
    """Find all polarishub tasks with leaderboards having at least 5 lines."""
    polarishub_dir = Path("biomlbench/tasks/polarishub")
    
    if not polarishub_dir.exists():
        print(f"Error: Directory {polarishub_dir} does not exist")
        return []
    
    populated_tasks = []
    
    # Find all leaderboard.csv files
    leaderboard_files = list(polarishub_dir.glob("*/leaderboard.csv"))
    
    print(f"Found {len(leaderboard_files)} leaderboard files in total")
    print("=" * 60)
    
    for leaderboard_path in sorted(leaderboard_files):
        task_name = leaderboard_path.parent.name
        line_count = count_non_empty_lines(leaderboard_path)
        
        print(f"Task: {task_name}")
        print(f"  Lines: {line_count}")
        
        if line_count >= 5:
            populated_tasks.append({
                'task_name': task_name,
                'line_count': line_count,
                'path': str(leaderboard_path)
            })
            print(f"  ✓ POPULATED (≥5 lines)")
        else:
            print(f"  ✗ Not populated (<5 lines)")
        print()
    
    return populated_tasks

def main():
    """Main function to run the script."""
    print("Finding PolarisHub tasks with populated leaderboards (≥5 lines)...")
    print()
    
    populated_tasks = find_populated_leaderboards()
    
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Tasks with populated leaderboards: {len(populated_tasks)}")
    print()
    
    if populated_tasks:
        print("Populated tasks:")
        for i, task in enumerate(populated_tasks, 1):
            print(f"{i:2d}. {task['task_name']} ({task['line_count']} lines)")
    else:
        print("No tasks found with leaderboards having ≥5 lines")

if __name__ == "__main__":
    main() 