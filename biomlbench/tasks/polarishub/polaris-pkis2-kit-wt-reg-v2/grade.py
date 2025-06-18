"""Auto-generated grading for polaris/pkis2-kit-wt-reg-v2."""
import pandas as pd
import polaris as po


def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """Grade using Polaris evaluation with main_metric."""
    benchmark = po.load_benchmark("polaris/pkis2-kit-wt-reg-v2")

    # Extract predictions (second column)
    predictions = submission[submission.columns[1]].tolist()

    # Use Polaris evaluation
    results = benchmark.evaluate(predictions)

    # Return main metric score
    main_metric = "mean_squared_error"

    # If main_metric is available in results, use it
    if hasattr(results, "results") and not results.results.empty:
        # Try to find main_metric in results
        if main_metric in results.results.columns:
            score = float(results.results[main_metric].iloc[0])
        else:
            # Fallback to first score column
            score_cols = [col for col in results.results.columns if col != "metric"]
            if score_cols:
                score = float(results.results[score_cols[0]].iloc[0])
            else:
                score = float(results.results.iloc[0, 0])
    else:
        # Fallback for older Polaris versions
        score = float(results.results.iloc[0]["score"])

    return score
