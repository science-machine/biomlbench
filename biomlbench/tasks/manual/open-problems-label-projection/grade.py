import numpy as np
import pandas as pd
from typing import List, Dict
from sklearn.metrics import f1_score
import sklearn.preprocessing

from biomlbench.grade_helpers import InvalidSubmissionError
from biomlbench.utils import get_logger

logger = get_logger(__name__)


def validate_submission_format(submission: pd.DataFrame, answers: pd.DataFrame) -> None:
    """
    Validate that the submission has the correct format for classification evaluation.
    
    Args:
        submission: Agent's predictions DataFrame
        answers: Ground truth answers DataFrame
        
    Raises:
        InvalidSubmissionError: If submission format is invalid
    """
    # Check required columns
    if 'cell_id' not in submission.columns:
        raise InvalidSubmissionError("Submission must contain 'cell_id' column")
    
    if 'label' not in submission.columns:
        raise InvalidSubmissionError("Submission must contain 'label' column")
    
    if len(submission) != len(answers):
        raise InvalidSubmissionError(
            f"Submission length ({len(submission)}) must match answers length ({len(answers)})"
        )
    
    # Check for duplicate cell IDs
    if submission['cell_id'].duplicated().any():
        raise InvalidSubmissionError("Submission contains duplicate cell IDs")
    
    # Check that all cell IDs in submission are present in answers
    submission_ids = set(submission['cell_id'])
    answer_ids = set(answers['cell_id'])
    
    missing_ids = answer_ids - submission_ids
    if missing_ids:
        raise InvalidSubmissionError(f"Submission missing {len(missing_ids)} cell IDs")
    
    extra_ids = submission_ids - answer_ids
    if extra_ids:
        raise InvalidSubmissionError(f"Submission contains {len(extra_ids)} extra cell IDs")
    
    # Check for null values
    if submission['label'].isnull().any():
        raise InvalidSubmissionError("Submission contains null labels")
    
    logger.info(f"Submission format validated: {len(submission)} cells")


def encode_labels(y_true: pd.Series, y_pred: pd.Series) -> tuple:
    """
    Encode categorical labels to integers using LabelEncoder.
    This ensures consistent encoding between true and predicted labels.
    
    Args:
        y_true: True labels (pandas Series)
        y_pred: Predicted labels (pandas Series)
        
    Returns:
        Tuple of (encoded_true, encoded_pred) as numpy arrays
    """
    # Get all unique labels from both true and predicted
    all_labels = list(set(y_true.unique()) | set(y_pred.unique()))
    
    # Create and fit encoder
    encoder = sklearn.preprocessing.LabelEncoder()
    encoder.fit(all_labels)
    
    # Transform labels
    y_true_encoded = encoder.transform(y_true)
    y_pred_encoded = encoder.transform(y_pred)
    
    logger.info(f"Encoded {len(encoder.classes_)} unique labels")
    
    return y_true_encoded, y_pred_encoded


def prepare_for_metric(submission: pd.DataFrame, answers: pd.DataFrame) -> tuple:
    """
    Prepare submission and answers data for F1 score calculation.
    
    Args:
        submission: Agent's predictions DataFrame
        answers: Ground truth answers DataFrame
        
    Returns:
        Tuple of (y_true, y_pred) as encoded arrays aligned by cell ID
    """
    # Sort both DataFrames by cell_id to ensure alignment
    submission_sorted = submission.sort_values('cell_id').reset_index(drop=True)
    answers_sorted = answers.sort_values('cell_id').reset_index(drop=True)
    
    # Verify alignment after sorting
    if not (submission_sorted['cell_id'] == answers_sorted['cell_id']).all():
        raise InvalidSubmissionError("Cell ID alignment failed after sorting")
    
    # Extract labels
    y_true = answers_sorted['label']
    y_pred = submission_sorted['label']
    
    # Encode labels using the same encoder
    y_true_encoded, y_pred_encoded = encode_labels(y_true, y_pred)
    
    logger.info(f"Prepared {len(y_true)} predictions for evaluation")
    
    return y_true_encoded, y_pred_encoded


def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """
    Grade a submission using F1-weighted score for cell type classification.
    
    This follows the OpenProblems benchmark approach:
    1. Validates submission format
    2. Encodes categorical labels consistently
    3. Computes F1-weighted score
    
    Args:
        submission: Agent's predictions DataFrame with 'cell_id' and 'label' columns
        answers: Ground truth DataFrame with 'cell_id' and 'label' columns
        
    Returns:
        F1-weighted score (float) - higher values are better (0 to 1)
    """
    try:
        # Validate submission format
        validate_submission_format(submission, answers)
        
        # Prepare data for evaluation (includes label encoding)
        y_true, y_pred = prepare_for_metric(submission, answers)
        
        # Calculate F1-weighted score
        score = f1_score(y_true, y_pred, average='weighted', zero_division=0)
        
        logger.info(f"F1-weighted score: {score:.4f}")
        
        # Additional metrics for logging (matching original benchmark)
        f1_macro = f1_score(y_true, y_pred, average='macro', zero_division=0)
        f1_micro = f1_score(y_true, y_pred, average='micro', zero_division=0)
        accuracy = np.mean(y_true == y_pred)
        
        logger.info(f"Additional metrics - F1-macro: {f1_macro:.4f}, "
                   f"F1-micro: {f1_micro:.4f}, Accuracy: {accuracy:.4f}")
        
        # Return the F1-weighted score (higher is better)
        return score
        
    except Exception as e:
        if isinstance(e, InvalidSubmissionError):
            raise
        else:
            raise InvalidSubmissionError(f"Error during evaluation: {str(e)}") from e 