#!/usr/bin/env python3
"""
Simple adapter to convert BioML-bench tasks to AutoBioML format.
Works with the actual container structure.
"""

import argparse
import sys
from pathlib import Path
import yaml


def read_instructions(data_dir: Path) -> str:
    """Read task instructions from mounted data."""
    # The main description file is description.md in the data directory
    desc_file = data_dir / "description.md"
    if desc_file.exists():
        return desc_file.read_text()
    else:
        raise FileNotFoundError(f"Description file not found in {data_dir}")


def detect_output_format(data_dir: Path) -> str:
    """Detect output format from sample submission file."""
    for pattern in ["sample_submission.*", "submission.*"]:
        sample_files = list(data_dir.glob(pattern))
        if sample_files:
            ext = sample_files[0].suffix.lower()
            if ext == ".csv":
                return "csv"
            elif ext == ".h5ad":
                return "h5ad"
            elif ext in [".arrow", ".feather"]:
                return "arrow"
            else:
                raise ValueError(f"Unsupported file format: {ext}")
        else:
            raise FileNotFoundError(f"Sample submission file not found in {data_dir}")


def create_task_yaml(task_id: str, data_dir: Path, work_dir: Path) -> Path:
    """Create a simple task.yaml file."""
    
    # Read instructions from mounted data
    description = read_instructions(data_dir)
    
    # Get list of data files
    files = [f.name for f in data_dir.glob("*") if f.is_file()]
    
    # Detect output format
    output_format = detect_output_format(data_dir)
    
    # Create simple task configuration
    # Use YAML literal scalar (|) to preserve formatting and avoid parsing issues
    task_config = {
        "name": task_id.replace("/", "_").replace("-", "_"),
        "display_name": f"BioML-bench Task: {task_id}",
        "version": "1.0",
        "task_description": description,  # This will be properly escaped by YAML
        "project_goal": "Complete the biomedical ML task",
        "available_data": {
            "agent_data": [
                {"path": f, "description": f"Data file: {f}"}
                for f in files
            ],
            "eval_data": []
        },
        "data_details": "See task description for details",
        "docker": {
            "gpu_required": False,
            "base_image": "millerh1/autobioml:latest"
        },
        # Add instructions directly but safely
        "full_instructions": f"{task_id}\n\n{description}\n\nData Files Available:\n" + 
                           "\n".join(f"- {f}" for f in files) + 
                           f"\n\nOutput Requirements:\nSave your predictions as `submission.{output_format}` in `/home/submission/`.\n\n" +
                           "Important Notes:\n" +
                           f"- Your submission MUST be saved as `/home/submission/submission.{output_format}`\n" +
                           "- Use the validation server at http://localhost:5000/validate to check your submission\n" +
                           "- Follow the task description above for specific requirements",
        "output_format": output_format
    }
    
    # Save task.yaml
    task_yaml_path = work_dir / "challenges" / "current" / "task.yaml"
    task_yaml_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(task_yaml_path, 'w') as f:
        yaml.dump(task_config, f, default_flow_style=False, allow_unicode=True)
    
    print(f"✓ Created task configuration")
    print(f"  - Output format: {output_format}")
    print(f"  - Data files: {len(files)}")
    
    return task_yaml_path


def setup_data_symlinks(data_dir: Path, work_dir: Path):
    """Create symlinks to data files."""
    for file in data_dir.glob("*"):
        if file.is_file():
            link_path = work_dir / "challenges" / "current" / file.name
            if not link_path.exists():
                link_path.symlink_to(file)


def main():
    parser = argparse.ArgumentParser(description="Simple BioML-bench to AutoBioML adapter")
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--submission-dir", required=True)
    
    args = parser.parse_args()
    
    # Use absolute paths based on container structure
    data_dir = Path("/home/data")  # Always mounted here
    work_dir = Path(args.work_dir)
    
    if not data_dir.exists():
        print(f"Error: Data directory {data_dir} does not exist")
        sys.exit(1)
    
    print(f"Container paths:")
    print(f"  Data dir: {data_dir}")
    print(f"  Work dir: {work_dir}")
    print(f"  Files in data: {list(data_dir.glob('*'))}")
    
    # Create task.yaml
    create_task_yaml(args.task_id, data_dir, work_dir)
    
    # Set up data symlinks
    setup_data_symlinks(data_dir, work_dir)
    
    # Ensure submission directory exists
    Path("/home/submission").mkdir(parents=True, exist_ok=True)
    
    print("Adapter completed")


if __name__ == "__main__":
    main() 