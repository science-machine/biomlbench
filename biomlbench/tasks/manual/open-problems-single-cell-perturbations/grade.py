import numpy as np
import pandas as pd
from typing import List

from biomlbench.grade_helpers import InvalidSubmissionError
from biomlbench.utils import get_logger

logger = get_logger(__name__)

# Expected number of genes (actual count from OpenProblems data)
EXPECTED_NUM_GENES = 5319


def get_gene_columns(df: pd.DataFrame) -> List[str]:
    """Extract gene column names from the dataframe."""
    # Exclude metadata columns to get just the gene columns
    metadata_cols = {'id', 'cell_type', 'sm_name', 'sm_lincs_id', 'SMILES', 'control'}
    gene_cols = [col for col in df.columns if col not in metadata_cols]
    return sorted(gene_cols)


def validate_submission_format(submission: pd.DataFrame, answers: pd.DataFrame) -> None:
    """
    Validate that the submission has the correct format for MRRMSE evaluation.
    
    Args:
        submission: Agent's predictions DataFrame
        answers: Ground truth answers DataFrame
        
    Raises:
        InvalidSubmissionError: If submission format is invalid
    """
    # Check required columns
    if 'id' not in submission.columns:
        raise InvalidSubmissionError("Submission must contain 'id' column")
    
    if len(submission) != len(answers):
        raise InvalidSubmissionError(
            f"Submission length ({len(submission)}) must match answers length ({len(answers)})"
        )
    
    # Extract gene columns from answers to validate against
    expected_gene_cols = get_gene_columns(answers)
    submission_gene_cols = get_gene_columns(submission)
    
    if len(expected_gene_cols) != EXPECTED_NUM_GENES:
        logger.warning(f"Expected {EXPECTED_NUM_GENES} genes but found {len(expected_gene_cols)} in answers")
    
    # Check that submission has all required gene columns
    missing_genes = set(expected_gene_cols) - set(submission_gene_cols)
    if missing_genes:
        raise InvalidSubmissionError(
            f"Submission is missing {len(missing_genes)} gene columns. "
            f"First 10 missing: {sorted(list(missing_genes))[:10]}"
        )
    
    extra_genes = set(submission_gene_cols) - set(expected_gene_cols)
    if extra_genes:
        raise InvalidSubmissionError(
            f"Submission has {len(extra_genes)} unexpected gene columns. "
            f"First 10 extra: {sorted(list(extra_genes))[:10]}"
        )
    
    # Check that all gene columns are numeric
    for gene in expected_gene_cols:
        if not pd.api.types.is_numeric_dtype(submission[gene]):
            raise InvalidSubmissionError(f"Gene column '{gene}' must be numeric")
    
    # Check for NaN values
    gene_data = submission[expected_gene_cols]
    if gene_data.isnull().any().any():
        null_counts = gene_data.isnull().sum()
        genes_with_nulls = null_counts[null_counts > 0]
        raise InvalidSubmissionError(
            f"Submission contains NaN values in {len(genes_with_nulls)} gene columns. "
            f"First 10: {list(genes_with_nulls.index[:10])}"
        )
    
    # Check that IDs match
    if not set(submission['id']) == set(answers['id']):
        missing_ids = set(answers['id']) - set(submission['id'])
        extra_ids = set(submission['id']) - set(answers['id'])
        error_msg = "Submission and answers must have the same IDs."
        if missing_ids:
            error_msg += f" Missing IDs: {sorted(list(missing_ids))[:10]}"
        if extra_ids:
            error_msg += f" Extra IDs: {sorted(list(extra_ids))[:10]}"
        raise InvalidSubmissionError(error_msg)


def mean_rowwise_rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """
    Calculate Mean Rowwise Root Mean Squared Error (MRRMSE).
    
    Formula: MRRMSE = (1/N_r) * sum(sqrt((1/N_c) * sum((y_true - y_pred)^2)))
    
    Where:
    - N_r is the number of rows (samples/cell_type-compound pairs)  
    - N_c is the number of columns (genes)
    - The inner sum is over columns (genes) for each row
    - The outer sum is over rows (samples)
    
    Args:
        y_true: Ground truth values, shape (n_samples, n_genes)
        y_pred: Predicted values, shape (n_samples, n_genes)
        
    Returns:
        MRRMSE score (lower is better)
    """
    if y_true.shape != y_pred.shape:
        raise ValueError(f"Shape mismatch: y_true {y_true.shape} vs y_pred {y_pred.shape}")
    
    # Calculate squared differences
    squared_diff = (y_true - y_pred) ** 2
    
    # Calculate mean squared error for each row (across genes)
    mse_per_row = np.mean(squared_diff, axis=1)
    
    # Calculate RMSE for each row
    rmse_per_row = np.sqrt(mse_per_row)
    
    # Calculate mean across all rows
    mrrmse = np.mean(rmse_per_row)
    
    return mrrmse


def prepare_for_metric(submission: pd.DataFrame, answers: pd.DataFrame) -> tuple:
    """
    Prepare submission and answers data for MRRMSE calculation.
    
    Args:
        submission: Agent's predictions DataFrame
        answers: Ground truth answers DataFrame
        
    Returns:
        Tuple of (y_true, y_pred) as numpy arrays
    """
    # Validate submission format
    validate_submission_format(submission, answers)
    
    # Get gene columns
    gene_cols = get_gene_columns(answers)
    
    # Sort both DataFrames by ID to ensure alignment
    submission_sorted = submission.sort_values('id').reset_index(drop=True)
    answers_sorted = answers.sort_values('id').reset_index(drop=True)
    
    # Verify ID alignment after sorting
    if not (submission_sorted['id'] == answers_sorted['id']).all():
        raise InvalidSubmissionError("ID alignment failed after sorting")
    
    # Extract gene expression values
    y_true = answers_sorted[gene_cols].values
    y_pred = submission_sorted[gene_cols].values
    
    logger.info(f"Prepared data for evaluation: {y_true.shape[0]} samples, {y_true.shape[1]} genes")
    
    return y_true, y_pred


def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """
    Grade a submission using Mean Rowwise Root Mean Squared Error.
    
    Args:
        submission: Agent's predictions DataFrame with 'id' and gene columns
        answers: Ground truth DataFrame with 'id' and gene columns
        
    Returns:
        MRRMSE score (float) - lower values are better
    """
    try:
        # Prepare data for evaluation
        y_true, y_pred = prepare_for_metric(submission, answers)

        # Clip values to -4 and 4
        y_pred = np.clip(y_pred, -4, 4)
        
        # Calculate MRRMSE
        score = mean_rowwise_rmse(y_true, y_pred)
        
        logger.info(f"MRRMSE score: {score:.6f}")
        
        # Return negative score since biomlbench expects higher scores to be better
        # but MRRMSE is a loss metric (lower is better)
        return score
        
    except Exception as e:
        if isinstance(e, InvalidSubmissionError):
            raise
        else:
            raise InvalidSubmissionError(f"Error during evaluation: {str(e)}") from e 