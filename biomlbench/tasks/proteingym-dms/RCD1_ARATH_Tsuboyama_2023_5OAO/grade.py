
"""Auto-generated grading for proteingym-dms/RCD1_ARATH_Tsuboyama_2023_5OAO."""
import pandas as pd
import polaris as po
from scipy.stats import spearmanr

def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """Grade using spearman correlation between true and predicted fitness scores."""

    # Extract predictions for each fold
    submission = submission[submission.iloc[:, 0] != "WT"]
    merged = pd.merge(submission, answers, on="id", how="inner")

    if len(merged) != len(submission):
        print(f"Warning: got a submission with {len(submission)} rows but {len(merged)} rows after merging with IDs with answers.")

    fold_random_5 = merged.iloc[:, 1].tolist()
    fold_modulo_5 = merged.iloc[:, 2].tolist()
    fold_contiguous_5 = merged.iloc[:, 3].tolist()
    answers = merged.iloc[:, 4].tolist()

    # Calculate Spearman correlation for each split
    correlation_fold_random_5, _ = spearmanr(fold_random_5, answers)
    correlation_fold_modulo_5, _ = spearmanr(fold_modulo_5, answers)
    correlation_fold_contiguous_5, _ = spearmanr(fold_contiguous_5, answers)

    # Return the average correlation across all splits
    return (correlation_fold_random_5 + correlation_fold_modulo_5 + correlation_fold_contiguous_5) / 3
