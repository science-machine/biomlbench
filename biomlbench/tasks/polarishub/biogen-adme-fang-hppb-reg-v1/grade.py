"""Auto-generated grading for biogen/adme-fang-hppb-reg-v1."""
import pandas as pd
import polaris as po

def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """Grade using Polaris evaluation with main_metric."""
    benchmark = po.load_benchmark("biogen/adme-fang-hppb-reg-v1")
    
    # Extract predictions (second column)
    predictions = submission[submission.columns[1]].tolist()
    
    # Use Polaris evaluation
    results = benchmark.evaluate(predictions)
    
    # Return main metric score
    main_metric = "pearsonr"
    
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
