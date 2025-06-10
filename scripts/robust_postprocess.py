#!/usr/bin/env python3
"""
Robust Post-Processing for BioML-bench Agent Outputs
Handles path corrections, environment issues, and validates outputs
"""

import os
import re
import sys
import shutil
import subprocess
from pathlib import Path
import pandas as pd
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

class RobustPostProcessor:
    def __init__(self, agent_dir, data_dir, output_dir):
        self.agent_dir = Path(agent_dir)
        self.data_dir = Path(data_dir)
        self.output_dir = Path(output_dir)
        self.corrections_applied = []
        
    def find_best_solution(self):
        """Find the best solution file with multiple fallback strategies"""
        search_paths = [
            self.agent_dir / "logs" / "best_solution.py",
            self.agent_dir / "best_solution.py",
            *list(self.agent_dir.glob("logs/**/best_solution.py")),
            *list(self.agent_dir.glob("logs/**/*.py")),
            *list(Path(".").glob("logs/**/*.py"))
        ]
        
        for path in search_paths:
            if path.exists() and path.stat().st_size > 0:
                logger.info(f"Found solution file: {path}")
                return path
                
        logger.warning("No solution file found")
        return None
    
    def apply_path_corrections(self, code_content):
        """Apply comprehensive path corrections to generated code"""
        original_content = code_content
        
        # Common path corrections
        corrections = [
            # Input directory fixes
            (r'"\.\/input\/', r'"./'),
            (r'"input\/', r'"./'),
            (r"'\.\/input\/", r"'./"),
            (r"'input\/", r"'./"),
            
            # Output directory fixes  
            (r'"\.\/working\/submission\.csv"', r'"submission.csv"'),
            (r'"working\/submission\.csv"', r'"submission.csv"'),
            (r'"\.\/output\/submission\.csv"', r'"submission.csv"'),
            (r'"\.\/submission\/submission\.csv"', r'"submission.csv"'),
            
            # to_csv path fixes
            (r'to_csv\("\.\/working\/', r'to_csv("'),
            (r'to_csv\("working\/', r'to_csv("'),
            (r'to_csv\("\.\/output\/', r'to_csv("'),
            
            # Directory creation removals (causes read-only filesystem errors)
            (r'os\.makedirs\("\.\/submission"[^)]*\)', r'# Directory creation removed'),
            (r'os\.makedirs\("\.\/output"[^)]*\)', r'# Directory creation removed'),
            (r'os\.makedirs\("\.\/working"[^)]*\)', r'# Directory creation removed'),
            (r'os\.makedirs\(.*submission.*\)', r'# Directory creation removed'),
            
            # Path object fixes
            (r'Path\("\.\/input"\)', r'Path(".")'),
            (r'Path\("input"\)', r'Path(".")'),
        ]
        
        for pattern, replacement in corrections:
            new_content = re.sub(pattern, replacement, code_content)
            if new_content != code_content:
                self.corrections_applied.append(f"Applied: {pattern} -> {replacement}")
                code_content = new_content
                
        return code_content
    
    def validate_data_files(self, work_dir):
        """Ensure required data files are available"""
        required_files = ['train.csv', 'test_features.csv']
        
        for filename in required_files:
            source = self.data_dir / filename
            target = work_dir / filename
            
            if source.exists():
                shutil.copy2(source, target)
                logger.info(f"Copied {filename} to working directory")
            else:
                logger.warning(f"Required file {filename} not found in {self.data_dir}")
                
        # Also copy sample submission if available
        sample_sub = self.data_dir / "sample_submission.csv"
        if sample_sub.exists():
            shutil.copy2(sample_sub, work_dir / "sample_submission.csv")
    
    def execute_solution(self, solution_file, work_dir):
        """Execute the solution in a controlled environment"""
        try:
            # Change to working directory
            original_cwd = os.getcwd()
            os.chdir(work_dir)
            
            # Execute the solution
            result = subprocess.run([
                sys.executable, str(solution_file)
            ], capture_output=True, text=True, timeout=600)  # 10 minute timeout
            
            logger.info(f"Solution execution completed with return code: {result.returncode}")
            
            if result.stdout:
                logger.info(f"STDOUT:\n{result.stdout}")
            if result.stderr:
                logger.warning(f"STDERR:\n{result.stderr}")
                
            return result.returncode == 0, result.stdout, result.stderr
            
        except subprocess.TimeoutExpired:
            logger.error("Solution execution timed out")
            return False, "", "Execution timed out"
        except Exception as e:
            logger.error(f"Solution execution failed: {e}")
            return False, "", str(e)
        finally:
            os.chdir(original_cwd)
    
    def validate_submission(self, submission_file):
        """Validate the generated submission file"""
        try:
            df = pd.read_csv(submission_file)
            
            # Check required columns
            required_columns = ['id', 'caco2_permeability']  # Task-specific
            missing_cols = [col for col in required_columns if col not in df.columns]
            if missing_cols:
                logger.error(f"Missing required columns: {missing_cols}")
                return False, f"Missing columns: {missing_cols}"
                
            # Check for realistic values (not all zeros)
            if 'caco2_permeability' in df.columns:
                pred_values = df['caco2_permeability'].values
                if len(set(pred_values)) == 1 and pred_values[0] == 0.0:
                    logger.warning("All predictions are zero - likely fallback data")
                    return False, "All predictions are zero"
                    
                # Check for reasonable range (log permeability typically -8 to -3)
                if pred_values.min() < -10 or pred_values.max() > 2:
                    logger.warning(f"Predictions outside reasonable range: {pred_values.min()} to {pred_values.max()}")
                
                # Check for NaN/infinite values
                if pd.isna(pred_values).any() or not pd.isfinite(pred_values).all():
                    logger.error("Predictions contain NaN or infinite values")
                    return False, "Invalid prediction values"
                    
            logger.info(f"Submission validation passed: {len(df)} predictions")
            return True, f"Valid submission with {len(df)} predictions"
            
        except Exception as e:
            logger.error(f"Submission validation failed: {e}")
            return False, str(e)
    
    def process(self):
        """Main processing pipeline"""
        logger.info("Starting robust post-processing...")
        
        # Find solution file
        solution_file = self.find_best_solution()
        if not solution_file:
            return False, "No solution file found"
            
        # Read and correct the solution code
        try:
            with open(solution_file, 'r') as f:
                original_code = f.read()
        except Exception as e:
            return False, f"Failed to read solution file: {e}"
            
        corrected_code = self.apply_path_corrections(original_code)
        
        # Create working directory
        work_dir = Path("/tmp") / f"biomlbench_work_{os.getpid()}"
        work_dir.mkdir(exist_ok=True)
        
        try:
            # Save corrected solution
            corrected_solution = work_dir / "corrected_solution.py"
            with open(corrected_solution, 'w') as f:
                f.write(corrected_code)
                
            # Copy data files
            self.validate_data_files(work_dir)
            
            # Execute solution
            success, stdout, stderr = self.execute_solution(corrected_solution, work_dir)
            
            # Check for submission file
            submission_files = list(work_dir.glob("submission.csv")) + list(work_dir.glob("**/submission.csv"))
            
            if submission_files:
                submission_file = submission_files[0]
                
                # Validate submission
                valid, message = self.validate_submission(submission_file)
                
                if valid:
                    # Copy to output location
                    output_submission = self.output_dir / "submission.csv"
                    shutil.copy2(submission_file, output_submission)
                    
                    logger.info(f"✅ Successfully processed submission: {message}")
                    logger.info(f"Applied {len(self.corrections_applied)} corrections")
                    
                    return True, f"Success: {message}"
                else:
                    logger.error(f"Invalid submission: {message}")
                    return False, f"Invalid submission: {message}"
            else:
                logger.error("No submission file generated")
                return False, "No submission file generated"
                
        finally:
            # Cleanup
            if work_dir.exists():
                shutil.rmtree(work_dir)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: robust_postprocess.py <agent_dir> <data_dir> <output_dir>")
        sys.exit(1)
        
    processor = RobustPostProcessor(sys.argv[1], sys.argv[2], sys.argv[3])
    success, message = processor.process()
    
    print(f"Result: {message}")
    sys.exit(0 if success else 1) 