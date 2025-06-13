"""
Utilities for biomedical task implementations.

This module contains helper functions commonly used across different biomedical tasks.
"""

import numpy as np
import pandas as pd

from biomlbench.grade_helpers import InvalidSubmissionError


def prepare_for_auroc_metric(
    submission: pd.DataFrame, 
    answers: pd.DataFrame, 
    id_col: str = "id", 
    target_col: str = "label"
) -> dict:
    """
    Prepare submission and answers for AUC-ROC calculation.
    
    Args:
        submission: DataFrame with predictions
        answers: DataFrame with ground truth
        id_col: Column name for IDs
        target_col: Column name for target values
        
    Returns:
        Dict with y_true and y_score arrays
    """
    # Merge on ID column
    merged = pd.merge(answers, submission, on=id_col, suffixes=('_true', '_pred'))
    
    y_true = merged[f"{target_col}_true"].values
    y_score = merged[f"{target_col}_pred"].values
    
    return {"y_true": y_true, "y_score": y_score}


def prepare_for_regression_metric(
    submission: pd.DataFrame, 
    answers: pd.DataFrame, 
    id_col: str = "id", 
    target_col: str = "label"
) -> dict:
    """
    Prepare submission and answers for regression metric calculation.
    
    Args:
        submission: DataFrame with predictions
        answers: DataFrame with ground truth
        id_col: Column name for IDs
        target_col: Column name for target values
        
    Returns:
        Dict with y_true and y_pred arrays
    """
    # Merge on ID column
    merged = pd.merge(answers, submission, on=id_col, suffixes=('_true', '_pred'))
    
    y_true = merged[f"{target_col}_true"].values
    y_pred = merged[f"{target_col}_pred"].values
    
    return {"y_true": y_true, "y_pred": y_pred}
