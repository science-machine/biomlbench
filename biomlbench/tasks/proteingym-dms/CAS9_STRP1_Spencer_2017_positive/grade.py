
"""Auto-generated grading for proteingym-dms/CAS9_STRP1_Spencer_2017_positive."""
import pandas as pd
import polaris as po
from scipy.stats import spearmanr

def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """Grade using spearman correlation between true and predicted fitness scores."""

    # Extract predictions (second column)
    predictions = submission[submission.columns[1]].tolist()
    answers = answers[answers.columns[1]].tolist()

    # Calculate Spearman correlation
    correlation, _ = spearmanr(predictions, answers)

    return correlation
