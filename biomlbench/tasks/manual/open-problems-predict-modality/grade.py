import numpy as np
import pandas as pd
from typing import List

from biomlbench.grade_helpers import InvalidSubmissionError
from biomlbench.utils import get_logger

logger = get_logger(__name__)


def get_protein_columns(df: pd.DataFrame) -> List[str]:
    """Extract protein column names from the dataframe."""
    # All columns except 'cell_id' are protein columns
    return [col for col in df.columns if col != 'cell_id']


def validate_submission_format(submission: pd.DataFrame, answers: pd.DataFrame) -> None:
    """
    Validate that the submission has the correct format for RMSE evaluation.
    
    Args:
        submission: Agent's predictions DataFrame
        answers: Ground truth answers DataFrame
        
    Raises:
        InvalidSubmissionError: If submission format is invalid
    """
    # Check required columns
    if 'cell_id' not in submission.columns:
        raise InvalidSubmissionError("Submission must contain 'cell_id' column")
    
    if 'cell_id' not in answers.columns:
        raise InvalidSubmissionError("Answers must contain 'cell_id' column")
    
    # Check submission length
    if len(submission) != len(answers):
        raise InvalidSubmissionError(
            f"Submission length ({len(submission)}) must match answers length ({len(answers)})"
        )
    
    # Check for duplicate cell IDs
    if submission['cell_id'].duplicated().any():
        raise InvalidSubmissionError("Submission contains duplicate cell IDs")
    
    if answers['cell_id'].duplicated().any():
        raise InvalidSubmissionError("Answers contains duplicate cell IDs")
    
    # Extract protein columns from answers
    expected_proteins = get_protein_columns(answers)
    submission_proteins = get_protein_columns(submission)
    
    # Check that submission has all required protein columns
    missing_proteins = set(expected_proteins) - set(submission_proteins)
    if missing_proteins:
        raise InvalidSubmissionError(
            f"Submission is missing {len(missing_proteins)} protein columns. "
            f"Missing: {sorted(list(missing_proteins))}"
        )
    
    extra_proteins = set(submission_proteins) - set(expected_proteins)
    if extra_proteins:
        raise InvalidSubmissionError(
            f"Submission has {len(extra_proteins)} unexpected protein columns. "
            f"Extra: {sorted(list(extra_proteins))}"
        )
    
    # Check that all protein columns are numeric
    for protein in expected_proteins:
        if not pd.api.types.is_numeric_dtype(submission[protein]):
            raise InvalidSubmissionError(f"Protein column '{protein}' must be numeric")
    
    # Check for NaN values
    protein_data = submission[expected_proteins]
    if protein_data.isnull().any().any():
        null_counts = protein_data.isnull().sum()
        proteins_with_nulls = null_counts[null_counts > 0]
        raise InvalidSubmissionError(
            f"Submission contains NaN values in {len(proteins_with_nulls)} protein columns. "
            f"Proteins with NaN: {list(proteins_with_nulls.index)}"
        )
    
    # Check that cell IDs match
    submission_cells = set(submission['cell_id'])
    answer_cells = set(answers['cell_id'])
    
    missing_cells = answer_cells - submission_cells
    if missing_cells:
        raise InvalidSubmissionError(
            f"Submission missing {len(missing_cells)} cell IDs. "
            f"First 10 missing: {sorted(list(missing_cells))[:10]}"
        )
    
    extra_cells = submission_cells - answer_cells
    if extra_cells:
        raise InvalidSubmissionError(
            f"Submission contains {len(extra_cells)} extra cell IDs. "
            f"First 10 extra: {sorted(list(extra_cells))[:10]}"
        )
    
    logger.info(f"Submission format validated: {len(submission)} cells, {len(expected_proteins)} proteins")


def rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """
    Calculate Root Mean Squared Error (RMSE).
    
    This matches the OpenProblems implementation:
    RMSE = sqrt(mean((y_true - y_pred)^2))
    
    The calculation is done across all values (all cells × all proteins).
    
    Args:
        y_true: Ground truth values, shape (n_cells, n_proteins)
        y_pred: Predicted values, shape (n_cells, n_proteins)
        
    Returns:
        RMSE score (lower is better)
    """
    if y_true.shape != y_pred.shape:
        raise ValueError(f"Shape mismatch: y_true {y_true.shape} vs y_pred {y_pred.shape}")
    
    # Calculate squared differences
    squared_diff = (y_true - y_pred) ** 2
    
    # Calculate mean across all values
    mse = np.mean(squared_diff)
    
    # Calculate RMSE
    rmse_value = np.sqrt(mse)
    
    return rmse_value


def prepare_for_metric(submission: pd.DataFrame, answers: pd.DataFrame) -> tuple:
    """
    Prepare submission and answers data for RMSE calculation.
    
    Args:
        submission: Agent's predictions DataFrame
        answers: Ground truth answers DataFrame
        
    Returns:
        Tuple of (y_true, y_pred) as numpy arrays
    """
    # Validate submission format
    validate_submission_format(submission, answers)
    
    # Get protein columns
    protein_cols = get_protein_columns(answers)
    
    # Sort both DataFrames by cell_id to ensure alignment
    submission_sorted = submission.sort_values('cell_id').reset_index(drop=True)
    answers_sorted = answers.sort_values('cell_id').reset_index(drop=True)
    
    # Verify cell_id alignment after sorting
    if not (submission_sorted['cell_id'] == answers_sorted['cell_id']).all():
        raise InvalidSubmissionError("Cell ID alignment failed after sorting")
    
    # Extract protein expression values
    y_true = answers_sorted[protein_cols].values
    y_pred = submission_sorted[protein_cols].values
    
    logger.info(f"Prepared data for evaluation: {y_true.shape[0]} cells, {y_true.shape[1]} proteins")
    
    return y_true, y_pred


def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """
    Grade a submission using Root Mean Squared Error (RMSE).
    
    Args:
        submission: Agent's predictions DataFrame with 'cell_id' and protein columns
        answers: Ground truth DataFrame with 'cell_id' and protein columns
        
    Returns:
        RMSE score (float) - lower values are better
        
    Note:
        The RMSE is calculated across all values (cells × proteins) to match
        the OpenProblems implementation.
    """
    try:
        # Prepare data for evaluation
        y_true, y_pred = prepare_for_metric(submission, answers)
        
        # Calculate RMSE
        score = rmse(y_true, y_pred)
        
        logger.info(f"RMSE score: {score:.6f}")
        
        # Return the RMSE score directly
        # Note: Lower is better for RMSE (it's a loss metric)
        return score
        
    except Exception as e:
        if isinstance(e, InvalidSubmissionError):
            raise
        else:
            raise InvalidSubmissionError(f"Error during evaluation: {str(e)}") from e 