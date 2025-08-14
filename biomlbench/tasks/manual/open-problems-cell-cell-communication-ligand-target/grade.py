import numpy as np
import pandas as pd
from biomlbench.grade_helpers import InvalidSubmissionError
from biomlbench.utils import get_logger

logger = get_logger(__name__)

def _sigmoid_transform(x):
    """Apply sigmoid transformation to normalize odds ratio."""
    return 1 - 1 / (1 + x / 2)

def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """
    Grade a submission using odds ratio metric.
    
    The odds ratio represents the ratio of true and false positives within a set
    of prioritized interactions (top ranked hits) versus the same ratio for the
    remainder of the interactions.
    
    Args:
        submission: DataFrame with columns [ligand, target, score]
        answers: DataFrame with columns [ligand, target, response, interaction_id]
    
    Returns:
        float: Odds ratio score after sigmoid transformation (0 to 1)
    """
    # Validate submission format
    required_columns = ["ligand", "target", "score"]
    if not all(col in submission.columns for col in required_columns):
        raise InvalidSubmissionError(
            f"Submission must contain columns: {required_columns}. "
            f"Found: {list(submission.columns)}"
        )
    
    # Check for missing values
    if submission.isnull().any().any():
        raise InvalidSubmissionError("Submission contains missing values")
    
    # Ensure score column is numeric
    try:
        submission["score"] = pd.to_numeric(submission["score"])
    except:
        raise InvalidSubmissionError("Score column must contain numeric values")
    
    # Check data types
    if not pd.api.types.is_numeric_dtype(submission["score"]):
        raise InvalidSubmissionError("Score column must be numeric")
    
    # Merge submission with ground truth (test set only)
    logger.info(f"Merging submission ({len(submission)} rows) with answers ({len(answers)} rows)")
    
    # Merge on ligand and target
    merged = pd.merge(
        answers[["ligand", "target", "response"]],
        submission[["ligand", "target", "score"]],
        on=["ligand", "target"],
        how="left"
    )
    
    # Check if all interactions have predictions
    missing_predictions = merged["score"].isnull().sum()
    if missing_predictions > 0:
        logger.warning(f"{missing_predictions} interactions missing predictions, assigning score of 0")
        merged["score"].fillna(0, inplace=True)
    
    # Sort by score in descending order
    merged = merged.sort_values("score", ascending=False)
    
    # Calculate odds ratio using top 5% of predictions
    top_prop = 0.05
    top_n = int(len(merged) * top_prop)
    
    if top_n < 1:
        top_n = 1  # At least one prediction in top set
    
    # Split into top and bottom sets
    top_set = merged.iloc[:top_n]
    bottom_set = merged.iloc[top_n:]
    
    # Calculate confusion matrix components
    tp = (top_set["response"] == 1).sum()  # True positives in top set
    fp = (top_set["response"] == 0).sum()  # False positives in top set
    fn = (bottom_set["response"] == 1).sum()  # False negatives in bottom set
    tn = (bottom_set["response"] == 0).sum()  # True negatives in bottom set
    
    logger.info(f"Confusion matrix - TP: {tp}, FP: {fp}, FN: {fn}, TN: {tn}")
    
    # Calculate odds ratio
    numerator = tp * tn
    denominator = fp * fn
    
    if denominator == 0:
        if numerator == 0:
            # Undefined case
            logger.warning("Odds ratio is undefined (0/0)")
            return 0.5  # Return neutral score
        else:
            # Perfect score
            logger.info("Perfect odds ratio (no false positives or false negatives)")
            odds_ratio = np.inf
    else:
        odds_ratio = numerator / denominator
    
    # Apply sigmoid transformation
    score = _sigmoid_transform(odds_ratio)
    
    logger.info(f"Odds ratio: {odds_ratio:.4f}, Transformed score: {score:.4f}")
    
    return score 