#!/usr/bin/env python3
"""
Adapter to convert BioML-bench tasks to AutoBioML format.
"""

import argparse
import json
import os
import shutil
import sys
from pathlib import Path
import yaml


def read_instructions(data_dir: Path) -> str:
    """Read task instructions from various possible locations."""
    # Check for instructions.txt in the parent directory
    instructions_file = data_dir.parent / "instructions.txt"
    if instructions_file.exists():
        return instructions_file.read_text()
    
    # Check for README or description files
    for filename in ["README.md", "readme.md", "description.txt", "task_description.txt"]:
        desc_file = data_dir / filename
        if desc_file.exists():
            return desc_file.read_text()
    
    return "Complete the biomedical machine learning task using the provided data."


def detect_task_type(data_dir: Path) -> dict:
    """Detect the task type and available data files."""
    files = list(data_dir.glob("*"))
    
    # Common patterns for different task types
    task_info = {
        "has_csv": any(f.suffix == ".csv" for f in files),
        "has_arrow": any(f.suffix == ".arrow" for f in files),
        "has_h5ad": any(f.suffix == ".h5ad" for f in files),
        "has_images": any(f.suffix in [".png", ".jpg", ".jpeg", ".tiff", ".dcm"] for f in files),
        "files": [f.name for f in files if f.is_file()]
    }
    
    # Detect specific task patterns
    if any("methylation" in f.name.lower() or "betas" in f.name.lower() for f in files):
        task_info["task_type"] = "epigenetic_clock"
    elif task_info["has_h5ad"]:
        task_info["task_type"] = "single_cell"
    elif task_info["has_images"]:
        task_info["task_type"] = "medical_imaging"
    else:
        task_info["task_type"] = "tabular_prediction"
    
    return task_info


def create_task_yaml(task_id: str, data_dir: Path, work_dir: Path) -> Path:
    """Create a task.yaml file compatible with AutoBioML."""
    
    instructions = read_instructions(data_dir)
    task_info = detect_task_type(data_dir)
    
    # Extract evaluation metrics from instructions if possible
    metrics = []
    if "pearson" in instructions.lower() or "correlation" in instructions.lower():
        metrics.append({
            "name": "pearson_correlation",
            "dataset": "test",
            "threshold": 0.7
        })
    elif "accuracy" in instructions.lower():
        metrics.append({
            "name": "accuracy",
            "dataset": "test",
            "threshold": 0.8
        })
    elif "auc" in instructions.lower() or "roc" in instructions.lower():
        metrics.append({
            "name": "auc_roc",
            "dataset": "test", 
            "threshold": 0.8
        })
    else:
        # Default metric
        metrics.append({
            "name": "custom_metric",
            "dataset": "test",
            "threshold": 0.5
        })
    
    # Create task configuration
    task_config = {
        "task_id": task_id,
        "display_name": f"BioML-bench Task: {task_id}",
        "version": "1.0.0",
        "description": instructions[:500] + "..." if len(instructions) > 500 else instructions,
        "data_dir": str(data_dir),
        "evaluation": {
            "metrics": metrics,
            "public_test_percentage": 0.5
        },
        "metadata": {
            "domain": "biomedical",
            "task_type": task_info["task_type"],
            "data_files": task_info["files"]
        }
    }
    
    # Save task.yaml
    task_yaml_path = work_dir / "challenges" / "current" / "task.yaml"
    task_yaml_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(task_yaml_path, 'w') as f:
        yaml.dump(task_config, f, default_flow_style=False)
    
    # Also create a more detailed instructions file
    instructions_path = work_dir / "challenges" / "current" / "instructions.md"
    with open(instructions_path, 'w') as f:
        f.write(f"# {task_id}\n\n")
        f.write(instructions)
        f.write("\n\n## Available Data Files\n\n")
        for file in task_info["files"]:
            f.write(f"- {file}\n")
    
    return task_yaml_path


def setup_data_symlinks(data_dir: Path, work_dir: Path):
    """Create symlinks to data files in the working directory."""
    # AutoBioML expects data files in the current working directory
    for file in data_dir.glob("*"):
        if file.is_file():
            link_path = work_dir / "challenges" / "current" / file.name
            if not link_path.exists():
                link_path.symlink_to(file)


def main():
    parser = argparse.ArgumentParser(description="Adapt BioML-bench task for AutoBioML")
    parser.add_argument("--task-id", required=True, help="Task ID")
    parser.add_argument("--data-dir", required=True, help="Data directory path")
    parser.add_argument("--work-dir", required=True, help="Working directory path")
    parser.add_argument("--submission-dir", required=True, help="Submission directory path")
    
    args = parser.parse_args()
    
    data_dir = Path(args.data_dir)
    work_dir = Path(args.work_dir)
    
    if not data_dir.exists():
        print(f"Error: Data directory {data_dir} does not exist")
        sys.exit(1)
    
    # Create task.yaml
    task_yaml_path = create_task_yaml(args.task_id, data_dir, work_dir)
    print(f"Created task configuration: {task_yaml_path}")
    
    # Set up data symlinks
    setup_data_symlinks(data_dir, work_dir)
    print(f"Set up data symlinks in {work_dir / 'challenges' / 'current'}")
    
    # Create submission directory
    submission_dir = Path(args.submission_dir)
    submission_dir.mkdir(parents=True, exist_ok=True)
    
    print("Task adaptation completed successfully")


if __name__ == "__main__":
    main() 