import numpy as np
import pandas as pd
from scipy.stats import kendalltau

from biomlbench.grade_helpers import InvalidSubmissionError
from biomlbench.utils import get_logger

logger = get_logger(__name__)


def validate_submission_format(submission: pd.DataFrame, answers: pd.DataFrame) -> None:
    """
    Validate that the submission has the correct format for Kendall correlation evaluation.
    
    Args:
        submission: Agent's predictions DataFrame
        answers: Ground truth answers DataFrame
        
    Raises:
        InvalidSubmissionError: If submission format is invalid
    """
    # Check required columns
    if 'gene_id' not in submission.columns:
        raise InvalidSubmissionError("Submission must contain 'gene_id' column")
    
    if 'spatial_score' not in submission.columns:
        raise InvalidSubmissionError("Submission must contain 'spatial_score' column")
    
    if len(submission) != len(answers):
        raise InvalidSubmissionError(
            f"Submission length ({len(submission)}) must match answers length ({len(answers)})"
        )
    
    # Check for duplicate gene IDs
    if submission['gene_id'].duplicated().any():
        raise InvalidSubmissionError("Submission contains duplicate gene IDs")
    
    # Check that all gene IDs in submission are present in answers
    submission_genes = set(submission['gene_id'])
    answer_genes = set(answers['gene_id'])
    
    missing_genes = answer_genes - submission_genes
    if missing_genes:
        raise InvalidSubmissionError(f"Submission missing {len(missing_genes)} gene IDs")
    
    extra_genes = submission_genes - answer_genes
    if extra_genes:
        raise InvalidSubmissionError(f"Submission contains {len(extra_genes)} extra gene IDs")
    
    # Validate spatial scores are numeric
    if not pd.api.types.is_numeric_dtype(submission['spatial_score']):
        raise InvalidSubmissionError("spatial_score column must contain numeric values")
    
    # Check for NaN values
    if submission['spatial_score'].isna().any():
        raise InvalidSubmissionError("spatial_score column contains NaN values")


def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """
    Grade a submission using Kendall's tau correlation coefficient.
    
    Args:
        submission: DataFrame with columns ['gene_id', 'spatial_score']
        answers: DataFrame with columns ['gene_id', 'true_spatial_var_score']
        
    Returns:
        Kendall's tau correlation coefficient (between -1 and 1, higher is better)
    """
    # Validate submission format
    validate_submission_format(submission, answers)
    
    # Merge submission with answers to ensure proper alignment
    merged = pd.merge(
        submission[['gene_id', 'spatial_score']],
        answers[['gene_id', 'true_spatial_var_score']],
        on='gene_id',
        how='inner'
    )
    
    if len(merged) != len(answers):
        raise InvalidSubmissionError("Failed to properly align submission with answers")
    
    # Extract predictions and true scores
    pred_scores = merged['spatial_score'].values
    true_scores = merged['true_spatial_var_score'].values
    
    # Calculate Kendall's tau correlation
    try:
        tau, p_value = kendalltau(pred_scores, true_scores)
        logger.info(f"Kendall's tau: {tau:.4f}, p-value: {p_value:.4e}")
        
        # Check if correlation is valid
        if np.isnan(tau):
            raise InvalidSubmissionError("Kendall's tau calculation resulted in NaN")
            
    except Exception as e:
        logger.error(f"Failed to compute Kendall's tau: {e}")
        raise InvalidSubmissionError(f"Cannot compute Kendall's tau: {e}")
    
    return tau


if __name__ == "__main__":
    # For testing
    import sys
    
    if len(sys.argv) != 3:
        print("Usage: python grade.py <submission_csv> <answers_csv>")
        sys.exit(1)
    
    submission = pd.read_csv(sys.argv[1])
    answers = pd.read_csv(sys.argv[2])
    
    try:
        score = grade(submission, answers)
        print(f"Kendall's tau: {score:.4f}")
    except InvalidSubmissionError as e:
        print(f"Invalid submission: {e}")
        sys.exit(1) 