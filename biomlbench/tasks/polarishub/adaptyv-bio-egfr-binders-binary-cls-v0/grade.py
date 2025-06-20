"""Auto-generated grading for adaptyv-bio/egfr-binders-binary-cls-v0."""
import pandas as pd
import polaris as po

def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """Grade using Polaris evaluation with main_metric."""
    benchmark = po.load_benchmark("adaptyv-bio/egfr-binders-binary-cls-v0")
    
    # Extract predictions (second column)
    predictions = submission[submission.columns[1]].tolist()
    
    # Use Polaris evaluation
    results = benchmark.evaluate(predictions)
    
    # Return main metric score
    main_metric = "balanced_accuracy"
    
    if hasattr(results, 'results') and not results.results.empty:
        # Try to find main_metric in results
        if "Metric" in results.results.columns:
            if results.results["Metric"].iloc[0] == main_metric:
                score = float(results.results["Score"].iloc[0])
            else:
                raise ValueError("This should never happen")
        else:
            raise ValueError("This should never happen")
    else:
        raise ValueError("This should never happen")
    
    return score
