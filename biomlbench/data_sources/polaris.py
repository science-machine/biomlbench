"""
Polaris data source implementation.

This module handles downloading datasets from Polaris Hub and provides
leaderboard information from the Polaris platform.
"""

import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd

from biomlbench.utils import get_logger
from biomlbench.data_sources.base import DataSource, DataSourceError
from .factory import register_data_source

logger = get_logger(__name__)


@register_data_source("polaris")
class PolarisDataSource(DataSource):
    """
    Data source for Polaris Hub benchmarks.
    
    Downloads benchmark data and provides leaderboard information
    from the Polaris platform using the polarishub conda environment.
    """
    
    def __init__(self):
        self.conda_env = "polarishub"  # Default environment name
    
    def validate_config(self, source_config: Dict[str, Any]) -> bool:
        """
        Validate Polaris source configuration.
        
        Args:
            source_config: Should contain 'benchmark_id' key
            
        Returns:
            True if valid
            
        Raises:
            DataSourceConfigError: If configuration is invalid
        """
        if 'benchmark_id' not in source_config:
            raise ValueError(
                "Polaris data source requires 'benchmark_id' in configuration",
                source_type="polaris"
            )
        
        benchmark_id = source_config['benchmark_id']
        if not isinstance(benchmark_id, str) or not benchmark_id.strip():
            raise ValueError(
                "Polaris 'benchmark_id' must be a non-empty string",
                source_type="polaris"
            )
        
        # Optional environment override
        if 'conda_env' in source_config:
            self.conda_env = source_config['conda_env']
        
        return True
    
    def _run_polaris_script(self, script_content: str) -> str:
        """
        Run a Python script in the polarishub conda environment.
        
        Args:
            script_content: Python script to execute
            
        Returns:
            Script output
            
        Raises:
            DataSourceError: If script execution fails
        """
        with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as f:
            f.write(script_content)
            script_path = f.name
        
        try:
            result = subprocess.run([
                'bash', '-c',
                f'source ~/miniconda3/etc/profile.d/conda.sh && '
                f'conda activate {self.conda_env} && '
                f'python {script_path}'
            ], capture_output=True, text=True, check=True)
            
            return result.stdout
            
        except subprocess.CalledProcessError as e:
            raise DataSourceError(
                f"Polaris script execution failed: {e.stderr}",
                source_type="polaris"
            ) from e
        finally:
            # Clean up temporary script
            Path(script_path).unlink()
    
    def download(self, source_config: Dict[str, Any], data_dir: Path) -> Optional[Path]:
        """
        Download benchmark data from Polaris Hub.
        
        Args:
            source_config: Must contain 'benchmark_id'
            data_dir: Directory to save data to
            
        Returns:
            Path to the data directory (Polaris doesn't use zip files)
            
        Raises:
            DataSourceError: If download fails
        """
        self.validate_config(source_config)
        
        benchmark_id = source_config['benchmark_id']
        data_dir.mkdir(parents=True, exist_ok=True)
        
        # Create script to download and save Polaris data
        download_script = f'''
import polaris as po
import pandas as pd
from pathlib import Path

try:
    # Load benchmark from Polaris Hub
    benchmark = po.load_benchmark('{benchmark_id}')
    train, _ = benchmark.get_train_test_split()
    # Grab the test set with the targets included
    test = benchmark._get_test_sets(hide_targets=False)['test']

    # Convert to dataframes
    df_train = train.as_dataframe()
    df_test = test.as_dataframe()

    # Save to data directory
    data_dir = Path('{data_dir}')
    df_train.to_csv(data_dir / 'polaris_train_data.csv', index=False)
    df_test.to_csv(data_dir / 'polaris_test_data.csv', index=False)

    print(f"SUCCESS: Downloaded train: {{len(df_train)}} samples")
    print(f"SUCCESS: Downloaded test: {{len(df_test)}} samples")
    print(f"SUCCESS: Saved to {{data_dir}}")
    
except Exception as e:
    print(f"ERROR: {{str(e)}}")
    raise
'''
        
        try:
            output = self._run_polaris_script(download_script)
            
            # Check if download was successful
            if "SUCCESS:" not in output:
                raise DataSourceError(
                    f"Polaris download failed: {output}",
                    source_type="polaris"
                )
            
            logger.info(f"Successfully downloaded Polaris benchmark '{benchmark_id}' to {data_dir}")
            return data_dir
            
        except Exception as e:
            if isinstance(e, DataSourceError):
                raise
            raise DataSourceError(
                f"Failed to download Polaris benchmark '{benchmark_id}': {e}",
                source_type="polaris"
            ) from e
    
    def get_leaderboard(self, source_config: Dict[str, Any]) -> pd.DataFrame:
        """
        Get leaderboard from Polaris Hub.
        
        Note: Polaris doesn't provide direct leaderboard access via API.
        This returns a minimal leaderboard structure.
        
        Args:
            source_config: Must contain 'benchmark_id'
            
        Returns:
            DataFrame with minimal leaderboard data
            
        Raises:
            DataSourceError: If leaderboard cannot be created
        """
        self.validate_config(source_config)
        
        benchmark_id = source_config['benchmark_id']
        
        try:
            # For now, create a basic leaderboard structure
            # In the future, this could scrape the Polaris Hub website
            # or use an API if one becomes available
            
            # Get benchmark info from Polaris
            info_script = f'''
import polaris as po

try:
    benchmark = po.load_benchmark('{benchmark_id}')
    metrics = getattr(benchmark, 'metrics', {{}})
    
    # Extract metric info
    if hasattr(metrics, '__iter__'):
        for metric in metrics:
            if hasattr(metric, 'name'):
                print(f"METRIC: {{metric.name}}")
                print(f"DIRECTION: {{getattr(metric, 'direction', 'unknown')}}")
                break
    
    print("SUCCESS: Retrieved benchmark info")
    
except Exception as e:
    print(f"ERROR: {{str(e)}}")
'''
            
            output = self._run_polaris_script(info_script)
            
            # Create a basic leaderboard structure
            # This should be enhanced with actual Polaris leaderboard data when available
            leaderboard_data = {
                'teamName': ['Polaris_Baseline', 'Random_Baseline'],
                'score': [0.5, 1.0],  # Placeholder scores
                'submissionDate': ['2024-01-01', '2024-01-01']
            }
            
            leaderboard_df = pd.DataFrame(leaderboard_data)
            
            logger.info(f"Created basic leaderboard structure for Polaris benchmark '{benchmark_id}'")
            return leaderboard_df
            
        except Exception as e:
            if isinstance(e, DataSourceError):
                raise
            raise DataSourceError(
                f"Failed to get leaderboard for Polaris benchmark '{benchmark_id}': {e}",
                source_type="polaris"
            ) from e 