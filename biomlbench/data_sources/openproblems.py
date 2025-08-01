"""Simple OpenProblems data source for PoC."""

import boto3
import yaml
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Any
import logging

from biomlbench.data_sources.base import DataSource
from biomlbench.data_sources.factory import register_data_source

logger = logging.getLogger(__name__)

@register_data_source("openproblems")
class OpenProblemsDataSource(DataSource):
    """Simple data source to extract OpenProblems results."""
    
    def __init__(self):
        self.bucket_name = "openproblems-data"
        self.s3_client = boto3.client('s3', region_name='us-east-1')
    
    def download(self, source_config: Dict[str, Any], data_dir: Path) -> Optional[Path]:
        """Not used for results-only tasks."""
        return None
    
    def get_leaderboard(self, source_config: Dict[str, Any]) -> pd.DataFrame:
        """Get leaderboard from OpenProblems S3."""
        task_name = source_config['task_name']  # Will raise KeyError if missing
        task_data = self.download_task_results(task_name)
        return self.create_simple_leaderboard(task_data)
    
    def download_task_results(self, task_name: str = "batch_integration") -> Dict:
        """Download results for batch_integration task."""
        # Hardcode the specific run we know exists
        run_prefix = f"resources/{task_name}/results/run_2024-06-28_13-20-27/"
        
        files_to_download = {
            'scores': 'score_uns.yaml',
            'metrics': 'metric_configs.yaml',
            'datasets': 'dataset_uns.yaml',
            'task_info': 'task_info.yaml'
        }
        
        task_data = {}
        
        for data_type, filename in files_to_download.items():
            s3_key = f"{run_prefix}{filename}"
            response = self.s3_client.get_object(Bucket=self.bucket_name, Key=s3_key)
            content = response['Body'].read().decode('utf-8')
            
            # Handle files with Groovy YAML tags by using unsafe loader or skipping problematic parts
            if filename in ['metric_configs.yaml', 'method_configs.yaml']:
                # These files contain Groovy-specific YAML that we can't safely parse
                # For PoC, we'll store as raw content and extract what we need later
                task_data[data_type] = {'_raw_content': content, '_parse_error': 'Groovy YAML tags not supported'}
                logger.info(f"Downloaded {filename} (stored as raw content due to Groovy YAML)")
            else:
                task_data[data_type] = yaml.safe_load(content)
                logger.info(f"Downloaded and parsed {filename}")
        
        return task_data
    
    def create_simple_leaderboard(self, task_data: Dict) -> pd.DataFrame:
        """Create a properly aggregated leaderboard with comparable scores."""
        scores = task_data['scores']  # Will raise KeyError if missing
        
        # Create one row per method-dataset-metric combination (full data)
        rows = []
        for entry in scores:
            method_id = entry['method_id']  # Will raise KeyError if missing
            dataset_id = entry['dataset_id']  # Will raise KeyError if missing
            metric_ids = entry['metric_ids']  # Will raise KeyError if missing
            metric_values = entry['metric_values']  # Will raise KeyError if missing
            date_created = entry['date_created']  # Will raise KeyError if missing
            
            # Create row for each metric
            for metric_id, metric_value in zip(metric_ids, metric_values):
                rows.append({
                    'method': method_id,
                    'dataset': dataset_id, 
                    'metric': metric_id,
                    'score': metric_value,
                    'date': date_created
                })
        
        full_df = pd.DataFrame(rows)
        
        # Create aggregated leaderboard: metric-balanced averaging
        # TODO: Revisit this -- we probably shouldn't just average scores across all datasets
        # 1. Average each method's score per metric across all datasets
        metric_averages = full_df.groupby(['method', 'metric'])['score'].mean().reset_index()
        
        # 2. Average across all metrics for each method (equal weight per metric)
        final_scores = metric_averages.groupby('method').agg({
            'score': 'mean',
            'metric': 'count'  # Count how many metrics each method has
        }).reset_index()
        
        # 3. Add date from the original data (use first occurrence)
        method_dates = full_df.groupby('method')['date'].first().reset_index()
        leaderboard = final_scores.merge(method_dates, on='method')
        
        # 4. Format for biomlbench leaderboard format
        leaderboard = leaderboard.rename(columns={
            'method': 'teamName',
            'score': 'score',
            'date': 'submissionDate'
        })
        
        # 5. Sort by score (descending) and return
        leaderboard = leaderboard.sort_values('score', ascending=False)
        
        return leaderboard[['teamName', 'score', 'submissionDate']]
    
    def get_detailed_results(self, source_config: Dict[str, Any]) -> pd.DataFrame:
        """Get detailed results showing all method-metric-dataset combinations."""
        task_name = source_config['task_name']  # Will raise KeyError if missing
        task_data = self.download_task_results(task_name)
        
        # Create one row per method-dataset-metric combination (full data)
        rows = []
        scores = task_data['scores']  # Will raise KeyError if missing
        
        for entry in scores:
            method_id = entry['method_id']  # Will raise KeyError if missing
            dataset_id = entry['dataset_id']  # Will raise KeyError if missing
            metric_ids = entry['metric_ids']  # Will raise KeyError if missing
            metric_values = entry['metric_values']  # Will raise KeyError if missing
            date_created = entry['date_created']  # Will raise KeyError if missing
            
            # Create row for each metric
            for metric_id, metric_value in zip(metric_ids, metric_values):
                rows.append({
                    'method': method_id,
                    'dataset': dataset_id, 
                    'metric': metric_id,
                    'score': metric_value,
                    'date': date_created
                })
        
        detailed_df = pd.DataFrame(rows)
        return detailed_df 