"""
Grading function for Caco-2 Wang permeability prediction task.

NOTE: This is a special case since Polaris benchmarks don't provide test labels.
For official evaluation, submissions should be sent to Polaris Hub.
This grader provides local evaluation for development/testing purposes.
"""

import pandas as pd

from biomlbench.tasks.utils import prepare_for_regression_metric


def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """
    Grade a Caco-2 permeability prediction submission using the Polaris evaluator.

    Args:
        submission: DataFrame with columns 'id' and 'caco2_permeability' (predictions)
        answers: DataFrame with columns 'id' and 'caco2_permeability' (true values - not used with Polaris)

    Returns:
        Mean Absolute Error from Polaris evaluation (lower is better)

    Note:
        This function uses the official Polaris benchmark for evaluation.
        The answers parameter is kept for compatibility but not used since
        Polaris handles the ground truth internally.
    """

    # Prepare data for regression evaluation to ensure proper ordering
    regression_inputs = prepare_for_regression_metric(
        submission=submission, answers=answers, id_col="id", target_col="caco2_permeability"
    )

    # Extract predictions in the correct order
    predictions = regression_inputs["y_pred"].tolist()

    # Use Polaris evaluator
    polaris_results = evaluate_with_polaris(predictions)

    return polaris_results["score"]


def evaluate_with_polaris(
    predictions: list[float], benchmark_name: str = "tdcommons/caco2-wang"
) -> dict:
    """
    Evaluate predictions using the official Polaris benchmark.

    Args:
        predictions: List of predicted Caco-2 permeability values
        benchmark_name: Name of the Polaris benchmark

    Returns:
        Dictionary containing Polaris evaluation results

    Note:
        This is the official way to evaluate on Polaris benchmarks.
    """
    try:
        import polaris as po
    except ImportError:
        raise ImportError(
            "Polaris is not installed. Please install it with: pip install polaris-hub"
        )

    # Load benchmark and evaluate
    benchmark = po.load_benchmark(benchmark_name)
    results = benchmark.evaluate(predictions)

    return {
        "score": float(results.results.Score.values[0]),
        "benchmark_name": benchmark_name,
        "num_predictions": len(predictions),
    }
