#!/usr/bin/env python3
"""
Convert AutoBioML output to BioML-bench submission format.
"""

import argparse
import shutil
import sys
from pathlib import Path
import pandas as pd
import numpy as np


def find_predictions_file(output_dir: Path) -> Path:
    """Find the predictions file in AutoBioML output."""
    # Common prediction file patterns
    patterns = [
        "predictions.arrow",
        "predictions.csv", 
        "submission.arrow",
        "submission.csv",
        "predictions_private.arrow",
        "final_predictions.arrow",
        "test_predictions.arrow"
    ]
    
    # Search in the output directory and subdirectories
    for pattern in patterns:
        # Direct match
        pred_file = output_dir / pattern
        if pred_file.exists():
            return pred_file
        
        # Search in subdirectories
        matches = list(output_dir.rglob(pattern))
        if matches:
            return matches[0]
    
    # Look for any arrow or csv file with predictions in the name
    arrow_files = list(output_dir.rglob("*prediction*.arrow"))
    if arrow_files:
        return arrow_files[0]
    
    csv_files = list(output_dir.rglob("*prediction*.csv"))
    if csv_files:
        return csv_files[0]
    
    return None


def find_h5ad_file(output_dir: Path) -> Path:
    """Find h5ad submission file for single-cell tasks."""
    patterns = [
        "submission.h5ad",
        "corrected.h5ad",
        "integrated.h5ad",
        "output.h5ad"
    ]
    
    for pattern in patterns:
        h5ad_file = output_dir / pattern
        if h5ad_file.exists():
            return h5ad_file
        
        matches = list(output_dir.rglob(pattern))
        if matches:
            return matches[0]
    
    # Look for any h5ad file
    h5ad_files = list(output_dir.rglob("*.h5ad"))
    if h5ad_files:
        return h5ad_files[0]
    
    return None


def convert_predictions(pred_file: Path, submission_dir: Path, task_id: str):
    """Convert predictions to the expected submission format."""
    
    # Determine file type and read
    if pred_file.suffix == ".arrow":
        df = pd.read_feather(pred_file)
    elif pred_file.suffix == ".csv":
        df = pd.read_csv(pred_file)
    else:
        raise ValueError(f"Unsupported file format: {pred_file.suffix}")
    
    # Handle different column naming conventions
    # AutoBioML might use different column names
    possible_pred_columns = ["prediction", "predictions", "pred", "y_pred", "target", "value"]
    pred_column = None
    
    for col in possible_pred_columns:
        if col in df.columns:
            pred_column = col
            break
    
    if pred_column is None and len(df.columns) == 1:
        # If there's only one column, assume it's the prediction
        pred_column = df.columns[0]
    
    if pred_column is None:
        # Try to find a numeric column that's not an ID
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        id_patterns = ["id", "index", "sample", "idx"]
        
        for col in numeric_cols:
            if not any(pattern in col.lower() for pattern in id_patterns):
                pred_column = col
                break
    
    if pred_column is None:
        raise ValueError(f"Could not identify prediction column in {df.columns.tolist()}")
    
    # Create submission dataframe
    submission_df = pd.DataFrame()
    
    # Check for ID column
    id_columns = ["id", "sample_id", "index", "sample_index", "idx"]
    id_column = None
    
    for col in id_columns:
        if col in df.columns:
            id_column = col
            break
    
    if id_column:
        submission_df["id"] = df[id_column]
    else:
        # Generate IDs if not present
        submission_df["id"] = range(len(df))
    
    # Add predictions
    submission_df["prediction"] = df[pred_column]
    
    # Save as CSV (BioML-bench standard format)
    submission_path = submission_dir / "submission.csv"
    submission_df.to_csv(submission_path, index=False)
    
    print(f"Created submission file: {submission_path}")
    print(f"Submission shape: {submission_df.shape}")
    print(f"Prediction column: {pred_column}")


def main():
    parser = argparse.ArgumentParser(description="Convert AutoBioML output to BioML-bench submission")
    parser.add_argument("--autobioml-output", required=True, help="AutoBioML output directory")
    parser.add_argument("--submission-dir", required=True, help="Submission directory")
    parser.add_argument("--task-id", required=True, help="Task ID")
    
    args = parser.parse_args()
    
    output_dir = Path(args.autobioml_output)
    submission_dir = Path(args.submission_dir)
    
    if not output_dir.exists():
        print(f"Error: Output directory {output_dir} does not exist")
        sys.exit(1)
    
    # Check for h5ad files first (single-cell tasks)
    h5ad_file = find_h5ad_file(output_dir)
    if h5ad_file:
        # Copy h5ad file directly
        submission_path = submission_dir / "submission.h5ad"
        shutil.copy2(h5ad_file, submission_path)
        print(f"Copied h5ad submission: {h5ad_file} -> {submission_path}")
        return
    
    # Look for predictions file
    pred_file = find_predictions_file(output_dir)
    
    if pred_file is None:
        print(f"Error: Could not find predictions file in {output_dir}")
        print("Searched for patterns: predictions.arrow, predictions.csv, submission.*, etc.")
        
        # List all files found
        all_files = list(output_dir.rglob("*"))
        print(f"\nFiles found in output directory:")
        for f in all_files:
            if f.is_file():
                print(f"  - {f.relative_to(output_dir)}")
        
        sys.exit(1)
    
    print(f"Found predictions file: {pred_file}")
    
    # Convert predictions
    try:
        convert_predictions(pred_file, submission_dir, args.task_id)
    except Exception as e:
        print(f"Error converting predictions: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 