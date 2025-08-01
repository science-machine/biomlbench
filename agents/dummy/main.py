"""
A dummy agent for single-cell batch integration tasks.

For CSV tasks: copies the sample_submission.csv or dataset.h5ad from the workspace/data directory.
For h5ad tasks: copies the input dataset.h5ad unchanged (no batch correction applied).
"""

import getpass
import os
import shutil
import sys
from pathlib import Path

# Get the current user's username
username = getpass.getuser()

# Check if the current user ID is 0 (root user ID on Unix-like systems)
if os.getuid() == 0:
    print(f"You are running this script as root. Your username is '{username}'.")
else:
    print(f"You do not have root access. Your username is {username}.")

print("The script is being run with the following python interpreter:")
print(sys.executable)

cwd = Path(__file__).parent
workspace_data_dir = cwd.parent / "data"
submission_dir = cwd.parent / "submission"

# Check what type of task this is based on available files
sample_submission_csv = workspace_data_dir / "sample_submission.csv"
dataset_h5ad = workspace_data_dir / "dataset.h5ad"

if dataset_h5ad.exists():
    # H5AD task (e.g., OpenProblems batch integration)
    print(f"Detected h5ad task - copying input dataset unchanged...")
    print(f"Input: {dataset_h5ad}")
    
    # For batch integration, the dummy agent should return the input unchanged
    # (i.e., no batch correction applied)
    output_path = submission_dir / "submission.h5ad"
    shutil.copy(dataset_h5ad, output_path)
    
    print(f"✅ Dummy submission created: {output_path}")
    print("Note: This submission applies NO batch correction (dummy baseline)")
    
elif sample_submission_csv.exists():
    # CSV task (traditional biomlbench tasks)
    print(f"Detected CSV task - copying sample submission...")
    
    output_path = submission_dir / "submission.csv"
    shutil.copy(sample_submission_csv, output_path)
    
    print(f"✅ Sample submission copied to {output_path}")
    
else:
    print("❌ Error: No sample_submission.csv or dataset.h5ad found!")
    print(f"Checked paths:")
    print(f"  - {sample_submission_csv}")
    print(f"  - {dataset_h5ad}")
    sys.exit(1)
