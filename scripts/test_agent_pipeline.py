#!/usr/bin/env python3
"""
BioML-bench Agent Pipeline Testing Framework
Comprehensive testing of agents before deployment
"""

import os
import sys
import json
import time
import shutil
import tempfile
import subprocess
from pathlib import Path
import pandas as pd
import logging
from typing import Dict, List, Tuple, Optional

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class AgentPipelineTester:
    def __init__(self, agent_name: str, test_timeout: int = 900):
        self.agent_name = agent_name
        self.test_timeout = test_timeout
        self.test_results = {}
        
    def create_test_data(self, test_dir: Path) -> bool:
        """Create minimal test dataset for validation"""
        try:
            # Create minimal but realistic test data
            train_data = {
                'smiles': [
                    'CCO',  # Ethanol
                    'CC(C)O',  # Isopropanol
                    'CCCCCCCCCCCCCCCCCC(=O)O',  # Stearic acid
                    'c1ccc(cc1)O',  # Phenol
                    'C1=CC=C(C=C1)C(=O)O'  # Benzoic acid
                ],
                'caco2_permeability': [-4.5, -4.2, -6.1, -4.8, -5.2]
            }
            
            test_data = {
                'id': [0, 1, 2],
                'smiles': ['CCCCO', 'c1ccc(cc1)CC', 'CC(C)(C)O']
            }
            
            sample_submission = {
                'id': [0, 1, 2],
                'caco2_permeability': [0.0, 0.0, 0.0]  # Template
            }
            
            # Save test files
            pd.DataFrame(train_data).to_csv(test_dir / 'train.csv', index=False)
            pd.DataFrame(test_data).to_csv(test_dir / 'test_features.csv', index=False)
            pd.DataFrame(sample_submission).to_csv(test_dir / 'sample_submission.csv', index=False)
            
            # Create answers for grading
            answers = {
                'id': [0, 1, 2],
                'caco2_permeability': [-4.3, -4.6, -4.1]  # Mock realistic answers
            }
            pd.DataFrame(answers).to_csv(test_dir / 'answers.csv', index=False)
            
            logger.info("✅ Test data created successfully")
            return True
            
        except Exception as e:
            logger.error(f"❌ Failed to create test data: {e}")
            return False
    
    def test_environment_compatibility(self) -> Tuple[bool, str]:
        """Test if agent environment is compatible with our requirements"""
        try:
            # Test basic imports
            test_script = """
import pandas as pd
import numpy as np
import sklearn
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

# Test basic molecular operations
mol = Chem.MolFromSmiles('CCO')
assert mol is not None
fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
assert len(fp) == 2048
mw = Descriptors.MolWt(mol)
assert mw > 0

print("Environment compatibility: PASS")
"""
            
            result = subprocess.run([
                'docker', 'run', '--rm', 
                f'{self.agent_name}:latest',
                'python', '-c', test_script
            ], capture_output=True, text=True, timeout=60)
            
            if result.returncode == 0:
                return True, "Environment compatibility test passed"
            else:
                return False, f"Environment test failed: {result.stderr}"
                
        except Exception as e:
            return False, f"Environment test error: {e}"
    
    def test_agent_execution(self, test_data_dir: Path) -> Tuple[bool, str, Dict]:
        """Test actual agent execution with monitoring"""
        try:
            # Create a test task set
            test_taskset = test_data_dir / "test_taskset.txt"
            with open(test_taskset, 'w') as f:
                f.write("caco2-wang\n")
            
            # Mount test data and run agent
            start_time = time.time()
            
            result = subprocess.run([
                'biomlbench', 'run-agent',
                '--agent', f'{self.agent_name}/dev',
                '--task-list', str(test_taskset)
            ], capture_output=True, text=True, timeout=self.test_timeout)
            
            execution_time = time.time() - start_time
            
            # Parse execution results
            execution_info = {
                'execution_time': execution_time,
                'return_code': result.returncode,
                'stdout_length': len(result.stdout),
                'stderr_length': len(result.stderr)
            }
            
            if result.returncode == 0:
                return True, f"Agent executed successfully in {execution_time:.1f}s", execution_info
            else:
                error_msg = result.stderr[-500:] if result.stderr else "Unknown error"
                return False, f"Agent execution failed: {error_msg}", execution_info
                
        except subprocess.TimeoutExpired:
            return False, f"Agent execution timed out after {self.test_timeout}s", {}
        except Exception as e:
            return False, f"Agent execution error: {e}", {}
    
    def test_output_validation(self, run_dir: Path) -> Tuple[bool, str, Dict]:
        """Validate agent outputs comprehensively"""
        try:
            validation_info = {}
            
            # Find submission file
            submission_files = list(run_dir.glob("**/submission.csv"))
            if not submission_files:
                return False, "No submission file found", validation_info
                
            submission_file = submission_files[0]
            validation_info['submission_file'] = str(submission_file)
            
            # Load and validate submission
            df = pd.read_csv(submission_file)
            validation_info['num_predictions'] = len(df)
            
            # Check structure
            required_cols = ['id', 'caco2_permeability']
            missing_cols = [col for col in required_cols if col not in df.columns]
            if missing_cols:
                return False, f"Missing columns: {missing_cols}", validation_info
                
            # Check prediction quality
            pred_values = df['caco2_permeability'].values
            validation_info['unique_values'] = len(set(pred_values))
            validation_info['all_zeros'] = len(set(pred_values)) == 1 and pred_values[0] == 0.0
            validation_info['value_range'] = [float(pred_values.min()), float(pred_values.max())]
            validation_info['has_nan'] = pd.isna(pred_values).any()
            validation_info['has_infinite'] = not pd.isfinite(pred_values).all()
            
            # Quality checks
            if validation_info['all_zeros']:
                return False, "All predictions are zero (likely fallback)", validation_info
                
            if validation_info['has_nan'] or validation_info['has_infinite']:
                return False, "Predictions contain invalid values", validation_info
                
            if validation_info['unique_values'] < 2:
                return False, "All predictions are identical", validation_info
                
            # Check for reasonable biological values
            min_val, max_val = validation_info['value_range']
            if min_val < -15 or max_val > 5:
                return False, f"Predictions outside biological range: {min_val:.2f} to {max_val:.2f}", validation_info
            
            # Check code quality
            code_files = list(run_dir.glob("**/best_solution.py"))
            if code_files:
                with open(code_files[0], 'r') as f:
                    code_content = f.read()
                    validation_info['uses_rdkit'] = 'rdkit' in code_content.lower()
                    validation_info['uses_ml'] = any(term in code_content.lower() 
                                                   for term in ['randomforest', 'xgboost', 'model.fit'])
                    validation_info['code_length'] = len(code_content)
            
            return True, f"Output validation passed: {validation_info['num_predictions']} realistic predictions", validation_info
            
        except Exception as e:
            return False, f"Output validation error: {e}", validation_info
    
    def test_performance_benchmarks(self, run_dir: Path) -> Tuple[bool, str, Dict]:
        """Test performance against expected benchmarks"""
        try:
            perf_info = {}
            
            # Check log files for performance metrics
            log_files = list(run_dir.glob("**/run.log"))
            if not log_files:
                return False, "No log file found", perf_info
                
            with open(log_files[0], 'r') as f:
                log_content = f.read()
                
            # Extract MAE/performance metrics
            import re
            mae_matches = re.findall(r'MAE[:\s]*([0-9.]+)', log_content)
            if mae_matches:
                mae_values = [float(m) for m in mae_matches]
                perf_info['mae_values'] = mae_values
                perf_info['best_mae'] = min(mae_values)
                
                # Performance benchmarks for Caco-2 permeability
                if perf_info['best_mae'] < 0.5:  # Reasonable performance
                    return True, f"Performance benchmark passed: MAE {perf_info['best_mae']:.4f}", perf_info
                else:
                    return False, f"Performance below benchmark: MAE {perf_info['best_mae']:.4f}", perf_info
            
            # Check execution time
            time_matches = re.findall(r'(\d+\.\d+) seconds', log_content)
            if time_matches:
                execution_time = float(time_matches[-1])
                perf_info['execution_time'] = execution_time
                
                if execution_time > 1800:  # 30 minutes max
                    return False, f"Execution too slow: {execution_time:.1f}s", perf_info
            
            return True, "Basic performance checks passed", perf_info
            
        except Exception as e:
            return False, f"Performance test error: {e}", perf_info
    
    def run_comprehensive_test(self) -> Dict:
        """Run all tests and return comprehensive results"""
        logger.info(f"🧪 Starting comprehensive test for agent: {self.agent_name}")
        
        overall_start = time.time()
        test_results = {
            'agent_name': self.agent_name,
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'tests': {},
            'overall_success': False
        }
        
        with tempfile.TemporaryDirectory() as temp_dir:
            test_dir = Path(temp_dir)
            
            # Test 1: Environment Compatibility
            logger.info("🔬 Testing environment compatibility...")
            success, message = self.test_environment_compatibility()
            test_results['tests']['environment'] = {
                'success': success,
                'message': message
            }
            
            if not success:
                logger.error(f"❌ Environment test failed: {message}")
                return test_results
                
            # Test 2: Create Test Data
            logger.info("📊 Creating test data...")
            if not self.create_test_data(test_dir):
                test_results['tests']['data_creation'] = {
                    'success': False,
                    'message': "Failed to create test data"
                }
                return test_results
                
            # Test 3: Agent Execution
            logger.info("🤖 Testing agent execution...")
            success, message, exec_info = self.test_agent_execution(test_dir)
            test_results['tests']['execution'] = {
                'success': success,
                'message': message,
                'details': exec_info
            }
            
            if not success:
                logger.error(f"❌ Agent execution failed: {message}")
                return test_results
                
            # Find the latest run directory
            run_dirs = list(Path("runs").glob("*_run-group_*"))
            if not run_dirs:
                test_results['tests']['no_run_dir'] = {
                    'success': False,
                    'message': "No run directory found"
                }
                return test_results
                
            latest_run = max(run_dirs, key=lambda x: x.stat().st_mtime)
            task_dirs = list(latest_run.glob("caco2-wang_*"))
            if not task_dirs:
                test_results['tests']['no_task_dir'] = {
                    'success': False,
                    'message': "No task directory found"
                }
                return test_results
                
            run_dir = task_dirs[0]
            
            # Test 4: Output Validation
            logger.info("📋 Validating outputs...")
            success, message, validation_info = self.test_output_validation(run_dir)
            test_results['tests']['output_validation'] = {
                'success': success,
                'message': message,
                'details': validation_info
            }
            
            # Test 5: Performance Benchmarks
            logger.info("⚡ Testing performance benchmarks...")
            success, message, perf_info = self.test_performance_benchmarks(run_dir)
            test_results['tests']['performance'] = {
                'success': success,
                'message': message,
                'details': perf_info
            }
            
        # Overall results
        all_tests_passed = all(
            test.get('success', False) 
            for test in test_results['tests'].values()
        )
        
        test_results['overall_success'] = all_tests_passed
        test_results['total_time'] = time.time() - overall_start
        
        if all_tests_passed:
            logger.info(f"🎉 All tests passed for {self.agent_name} in {test_results['total_time']:.1f}s")
        else:
            failed_tests = [name for name, test in test_results['tests'].items() if not test.get('success', False)]
            logger.error(f"❌ Tests failed for {self.agent_name}: {failed_tests}")
            
        return test_results

def main():
    if len(sys.argv) != 2:
        print("Usage: test_agent_pipeline.py <agent_name>")
        sys.exit(1)
        
    agent_name = sys.argv[1]
    tester = AgentPipelineTester(agent_name)
    results = tester.run_comprehensive_test()
    
    # Save results
    results_file = f"test_results_{agent_name}_{int(time.time())}.json"
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
        
    logger.info(f"📊 Test results saved to: {results_file}")
    
    # Exit with appropriate code
    sys.exit(0 if results['overall_success'] else 1)

if __name__ == "__main__":
    main() 