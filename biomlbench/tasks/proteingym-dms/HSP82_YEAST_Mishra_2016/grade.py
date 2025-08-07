
"""Auto-generated grading for proteingym-dms/HSP82_YEAST_Mishra_2016."""
import pandas as pd
import polaris as po
from scipy.stats import spearmanr

def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """Grade using spearman correlation between true and predicted fitness scores."""

    # Extract predictions for each fold
    fold_random_5 = submission[submission.columns[1]].tolist()
    fold_modulo_5 = submission[submission.columns[2]].tolist()
    fold_contiguous_5 = submission[submission.columns[3]].tolist()
    answers_fold_random_5 = answers[answers.columns[1]].tolist()
    answers_fold_modulo_5 = answers[answers.columns[2]].tolist()
    answers_fold_contiguous_5 = answers[answers.columns[3]].tolist()

    # Calculate Spearman correlation for each split
    correlation_fold_random_5, _ = spearmanr(fold_random_5, answers_fold_random_5)
    correlation_fold_modulo_5, _ = spearmanr(fold_modulo_5, answers_fold_modulo_5)
    correlation_fold_contiguous_5, _ = spearmanr(fold_contiguous_5, answers_fold_contiguous_5)

    # Return the average correlation across all splits
    return (correlation_fold_random_5 + correlation_fold_modulo_5 + correlation_fold_contiguous_5) / 3
