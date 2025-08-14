"""Auto-generated grading for tdcommons/cyp2d6-substrate-carbonmangels."""
import pandas as pd
import polaris as po

def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """Grade using Polaris evaluation with main_metric."""
    benchmark = po.load_benchmark("tdcommons/cyp2d6-substrate-carbonmangels")
    
    # Extract predictions (second column)
    predictions = submission[submission.columns[1]].tolist()
    
    # Use Polaris evaluation
    try:
        results = benchmark.evaluate(predictions)
    except ValueError as e:
        if "Metric requires `y_prob` input" in str(e):
            results = benchmark.evaluate(y_prob=predictions)
        else:
            raise e
    
    # Return main metric score
    main_metric = "pr_auc"
    
    if hasattr(results, 'results') and not results.results.empty:
        # Search for the main metric in the results
        if "Metric" in results.results.columns:
            # Find the row that matches the main metric
            metric_mask = results.results["Metric"] == main_metric
            matching_rows = results.results[metric_mask]
            
            if matching_rows.empty:
                raise ValueError(
                    f"Main metric 'pr_auc' not found in evaluation results. "
                )
            elif len(matching_rows) > 1:
                raise ValueError(
                    f"Multiple rows found for main metric 'pr_auc'. "
                )
            else:
                # Extract score from the single matching row
                score = float(matching_rows["Score"].iloc[0])
        else:
            raise ValueError(
                f"'Metric' column not found in evaluation results. "
            )
    else:
        raise ValueError(
            "Evaluation results are empty or missing. "
        )
    
    return score
