"""
Grading function for Caco-2 Wang permeability prediction task.

NOTE: This is a special case since Polaris benchmarks don't provide test labels.
For official evaluation, submissions should be sent to Polaris Hub.
This grader provides local evaluation for development/testing purposes.
"""

import pandas as pd
from sklearn.metrics import mean_absolute_error

from biomlbench.tasks.utils import prepare_for_regression_metric


def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """
    Grade a Caco-2 permeability prediction submission using Mean Absolute Error.
    
    Args:
        submission: DataFrame with columns 'id' and 'caco2_permeability' (predictions)
        answers: DataFrame with columns 'id' and 'caco2_permeability' (true values)
        
    Returns:
        Mean Absolute Error (lower is better)
        
    Note:
        Since Polaris benchmarks don't provide test labels, this function
        will only work with synthetic/validation data. For official evaluation,
        use the Polaris Hub submission system.
    """
    
    # Prepare data for regression evaluation
    regression_inputs = prepare_for_regression_metric(
        submission=submission, 
        answers=answers, 
        id_col="id", 
        target_col="caco2_permeability"
    )
    
    # Calculate Mean Absolute Error
    mae = mean_absolute_error(
        y_true=regression_inputs["y_true"], 
        y_pred=regression_inputs["y_pred"]
    )
    
    return mae


def evaluate_with_polaris(predictions: list[float], benchmark_name: str = "tdcommons/caco2-wang") -> dict:
    """
    Evaluate predictions using the official Polaris benchmark.
    
    This function requires the polarishub conda environment and Polaris authentication.
    
    Args:
        predictions: List of predicted Caco-2 permeability values
        benchmark_name: Name of the Polaris benchmark
        
    Returns:
        Dictionary containing Polaris evaluation results
        
    Note:
        This is the official way to evaluate on Polaris benchmarks.
        The local grade() function is only for development purposes.
    """
    import subprocess
    import json
    import tempfile
    
    # Create temporary file with predictions
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        json.dump(predictions, f)
        pred_file = f.name
    
    # Create evaluation script
    eval_script = f"""
import polaris as po
import json

# Load predictions
with open('{pred_file}', 'r') as f:
    predictions = json.load(f)

# Load benchmark and evaluate
benchmark = po.load_benchmark('{benchmark_name}')
results = benchmark.evaluate(predictions)

# Print results
print("POLARIS_RESULTS_START")
print(json.dumps({{
    'score': float(results.scores['mean_absolute_error']),
    'benchmark_name': '{benchmark_name}',
    'num_predictions': len(predictions)
}}))
print("POLARIS_RESULTS_END")
"""
    
    script_path = 'temp_polaris_eval.py'
    with open(script_path, 'w') as f:
        f.write(eval_script)
    
    try:
        # Run evaluation in polarishub environment
        result = subprocess.run([
            'bash', '-c',
            f'source ~/miniconda3/etc/profile.d/conda.sh && conda activate polarishub && python {script_path}'
        ], capture_output=True, text=True, check=True)
        
        # Extract results from output
        output_lines = result.stdout.split('\n')
        start_idx = None
        end_idx = None
        
        for i, line in enumerate(output_lines):
            if "POLARIS_RESULTS_START" in line:
                start_idx = i + 1
            elif "POLARIS_RESULTS_END" in line:
                end_idx = i
                break
        
        if start_idx is not None and end_idx is not None:
            results_json = '\n'.join(output_lines[start_idx:end_idx])
            return json.loads(results_json)
        else:
            raise ValueError("Could not parse Polaris results from output")
            
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Polaris evaluation failed: {e.stderr}")
        
    finally:
        # Clean up temporary files
        import os
        os.remove(pred_file)
        os.remove(script_path) 