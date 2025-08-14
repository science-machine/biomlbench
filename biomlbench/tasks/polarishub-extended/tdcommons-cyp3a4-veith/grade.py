"""Auto-generated grading for tdcommons/cyp3a4-veith."""
import pandas as pd
import polaris as po

def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """Grade using Polaris evaluation with main_metric."""
    benchmark = po.load_benchmark("tdcommons/cyp3a4-veith")
    
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
