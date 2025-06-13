import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from biomlbench.tasks.utils import prepare_for_auroc_metric


def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """
    Grade a histopathologic cancer detection submission using AUC-ROC.
    
    Args:
        submission: DataFrame with columns 'id' and 'label' (predicted probabilities)
        answers: DataFrame with columns 'id' and 'label' (true binary labels)
        
    Returns:
        AUC-ROC score between 0 and 1
    """
    roc_auc_inputs = prepare_for_auroc_metric(
        submission=submission, answers=answers, id_col="id", target_col="label"
    )
    return roc_auc_score(y_true=roc_auc_inputs["y_true"], y_score=roc_auc_inputs["y_score"])
