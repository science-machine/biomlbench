
"""Auto-generated grading for proteingym-dms/RAD_ANTMA_Tsuboyama_2023_2CJJ_indels."""
import pandas as pd
import polaris as po
from scipy.stats import spearmanr

def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """Grade using spearman correlation between true and predicted fitness scores."""

    # Extract predictions for each fold
    submission = submission[submission.iloc[:, 0] != "WT"]
    merged = pd.merge(submission, answers, on="id", how="inner", suffixes=("_pred", "_true"))

    if len(merged) != len(submission):
        print(f"Warning: got a submission with {len(submission)} rows but {len(merged)} rows after merging with IDs with answers.")

    correlation, _ = spearmanr(merged["fitness_score_pred"], merged["fitness_score_true"])

    return correlation
