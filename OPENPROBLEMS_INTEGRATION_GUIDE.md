# OpenProblems Integration Guide for biomlbench

## ⚠️ CRITICAL REVIEW AND CORRECTIONS

After cross-referencing with OpenProblems documentation, GitHub repositories, and actual YAML data structures, several critical issues have been identified in the original approach. **This section provides the corrected implementation strategy.**

### **Major Issues Identified**

1. **Architecture Mismatch**: OpenProblems uses Viash components with AnnData files, not standalone Python scripts
2. **Data Format Incompatibility**: Missing integration with AnnData (.h5ad) files 
3. **Component Structure**: Ignoring OpenProblems' component-based architecture
4. **Metric Handling**: Incorrectly aggregating metrics instead of preserving individual scores
5. **Evaluation Approach**: Placeholder scripts that don't actually work with OpenProblems data

### **Corrected Implementation Approach**

#### **Option 1: Results-Only Integration (Recommended)**

Instead of trying to recreate OpenProblems tasks, **extract and present benchmark results only**:

```python
class OpenProblemsResultsSource:
    """Extract benchmark results from OpenProblems for display/analysis."""
    
    def get_leaderboard_data(self, task_name: str) -> pd.DataFrame:
        """Get formatted leaderboard preserving all metrics."""
        # Download results from S3
        scores = self.download_scores(task_name)
        
        # Create leaderboard preserving individual metrics
        rows = []
        for entry in scores:
            for metric_id, metric_value in zip(entry['metric_ids'], entry['metric_values']):
                rows.append({
                    'task': task_name,
                    'method': entry['method_id'],
                    'dataset': entry['dataset_id'],
                    'metric': metric_id,
                    'score': metric_value,
                    'date': entry['date_created']
                })
        
        return pd.DataFrame(rows)
    
    def generate_summary_report(self, task_name: str) -> Dict:
        """Generate summary report for a task."""
        return {
            'task_info': self.get_task_info(task_name),
            'leaderboard': self.get_leaderboard_data(task_name),
            'metrics_info': self.get_metrics_info(task_name),
            'datasets_info': self.get_datasets_info(task_name)
        }
```

#### **Option 2: Component Bridge Integration (Advanced)**

Create a bridge to execute actual OpenProblems components:

```python
class OpenProblemsComponentBridge:
    """Bridge to execute OpenProblems Viash components."""
    
    def setup_task_environment(self, task_name: str):
        """Set up OpenProblems task environment."""
        # Clone task repository
        # Install dependencies (Viash, Nextflow, etc.)
        # Download required datasets
        pass
    
    def run_method(self, task_name: str, method_name: str, input_data: str) -> str:
        """Execute an OpenProblems method component."""
        # Use Viash to run the component
        # Return path to output file
        pass
    
    def evaluate_results(self, task_name: str, predictions: str, solution: str) -> Dict:
        """Evaluate using OpenProblems metrics."""
        # Run OpenProblems metric components
        # Return metric scores
        pass
```

#### **Option 3: Hybrid Approach (Balanced)**

Combine results extraction with simplified task recreation:

```python
class OpenProblemsHybridIntegration:
    """Hybrid approach: Extract results + create reference tasks."""
    
    def create_reference_task(self, task_name: str):
        """Create a reference task with actual OpenProblems data."""
        
        # 1. Extract benchmark results for leaderboard
        results = self.extract_results(task_name)
        
        # 2. Create a prepare.py that points to OpenProblems datasets
        prepare_script = self.generate_data_pointer_script(task_name)
        
        # 3. Create a grade.py that explains OpenProblems metrics
        grade_script = self.generate_reference_grader(task_name)
        
        # 4. Include links to run actual OpenProblems evaluation
        instructions = self.generate_evaluation_instructions(task_name)
        
        return {
            'config.yaml': self.generate_config(task_name, results),
            'leaderboard.csv': self.format_leaderboard(results),
            'prepare.py': prepare_script,
            'grade.py': grade_script,
            'evaluation_guide.md': instructions,
            'description.md': self.generate_description(task_name)
        }
```

### **Recommended Corrected Implementation**

Based on the analysis, **Option 1 (Results-Only Integration)** is most appropriate because:

1. **Realistic Scope**: Focuses on what's actually achievable
2. **Value Preservation**: Maintains all OpenProblems metric information  
3. **Compatibility**: Works with existing biomlbench architecture
4. **Maintainability**: Doesn't require complex Viash/Nextflow integration

### **Updated Phase 1: Results Extraction Data Source**

```python
class OpenProblemsResultsSource(BaseDataSource):
    """Extract and format OpenProblems benchmark results."""
    
    def __init__(self):
        self.bucket_name = "openproblems-data"
        self.base_prefix = "resources"
        self.s3_client = boto3.client('s3', region_name='us-east-1')
    
    def get_task_results(self, task_name: str) -> Dict:
        """Get complete results data for a task."""
        # Download all YAML files
        run_prefix = self.get_latest_run(task_name)
        
        return {
            'scores': self.download_yaml(f"{run_prefix}score_uns.yaml"),
            'metrics': self.download_yaml(f"{run_prefix}metric_configs.yaml"),
            'datasets': self.download_yaml(f"{run_prefix}dataset_uns.yaml"),
            'task_info': self.download_yaml(f"{run_prefix}task_info.yaml"),
            'methods': self.download_yaml(f"{run_prefix}method_configs.yaml"),
            'state': self.download_yaml(f"{run_prefix}state.yaml")
        }
    
    def create_detailed_leaderboard(self, task_name: str) -> pd.DataFrame:
        """Create detailed leaderboard preserving all metrics."""
        results = self.get_task_results(task_name)
        
        # Build comprehensive leaderboard
        leaderboard_data = []
        for score_entry in results['scores']:
            method_id = score_entry.get('method_id', 'unknown')
            dataset_id = score_entry.get('dataset_id', 'unknown')
            date_created = score_entry.get('date_created', '2024-01-01')
            
            # Preserve each metric separately
            metric_ids = score_entry.get('metric_ids', [])
            metric_values = score_entry.get('metric_values', [])
            
            for metric_id, metric_value in zip(metric_ids, metric_values):
                leaderboard_data.append({
                    'task': task_name,
                    'method': method_id,
                    'dataset': dataset_id,
                    'metric': metric_id,
                    'score': metric_value,
                    'submission_date': date_created,
                    'normalization': score_entry.get('normalization_id', 'unknown')
                })
        
        return pd.DataFrame(leaderboard_data)
    
    def generate_task_summary(self, task_name: str) -> Dict:
        """Generate comprehensive task summary."""
        results = self.get_task_results(task_name)
        leaderboard = self.create_detailed_leaderboard(task_name)
        
        # Extract key statistics
        summary = {
            'task_name': task_name,
            'total_methods': len(leaderboard['method'].unique()),
            'total_datasets': len(leaderboard['dataset'].unique()),
            'total_metrics': len(leaderboard['metric'].unique()),
            'total_evaluations': len(leaderboard),
            'date_range': {
                'earliest': leaderboard['submission_date'].min(),
                'latest': leaderboard['submission_date'].max()
            },
            'top_methods_by_metric': self.get_top_methods_by_metric(leaderboard),
            'metric_definitions': self.extract_metric_definitions(results['metrics']),
            'dataset_descriptions': self.extract_dataset_info(results['datasets']),
            'task_description': results['task_info'],
            'evaluation_setup': results['state']
        }
        
        return summary
```

## Original Implementation (DEPRECATED)

**⚠️ The following implementation has been identified as problematic and should not be used as-is. It's preserved for reference only.**

## Overview

This document provides comprehensive step-by-step instructions for integrating OpenProblems benchmarks into the biomlbench framework. OpenProblems is a community-driven benchmarking platform for single-cell analysis methods that stores comprehensive benchmark results in AWS S3.

## Table of Contents

1. [Prerequisites and Setup](#prerequisites-and-setup)
2. [Understanding the Data Structures](#understanding-the-data-structures)
3. [Creating the OpenProblems Data Source](#creating-the-openproblems-data-source)
4. [Building the Ingestion Pipeline](#building-the-ingestion-pipeline)
5. [Data Format Conversion](#data-format-conversion)
6. [Task Generation](#task-generation)
7. [Testing and Validation](#testing-and-validation)
8. [Integration with biomlbench](#integration-with-biomlbench)
9. [Maintenance and Updates](#maintenance-and-updates)

## Prerequisites and Setup

### Required Dependencies

Add these to `pyproject.toml`:

```toml
[project.dependencies]
boto3 = ">=1.26.0"
s3fs = ">=2023.1.0" 
pyyaml = ">=6.0"
```

### AWS Configuration

OpenProblems data is stored in a public S3 bucket. Configure boto3:

```python
# No credentials needed for public bucket, but boto3 must be available
import boto3
s3 = boto3.client('s3', region_name='us-east-1')
```

## Understanding the Data Structures

### OpenProblems S3 Structure

```
s3://openproblems-data/resources/
├── batch_integration/
│   └── results/
│       └── run_2024-06-28_13-20-27/
│           ├── score_uns.yaml          # Leaderboard data
│           ├── metric_configs.yaml     # Metric definitions
│           ├── method_configs.yaml     # Method configurations
│           ├── dataset_uns.yaml        # Dataset metadata
│           ├── task_info.yaml         # Task description
│           └── state.yaml             # Run metadata
├── label_projection/
├── perturbation_prediction/
└── [other tasks]/
```

### Key Data File Formats

**score_uns.yaml**: Contains benchmark results
```yaml
- dataset_id: cellxgene_census/immune_cell_atlas
  date_created: 28-06-2024
  method_id: fastmnn_embedding
  metric_ids: [asw_label]
  metric_values: [0.5377879595923591]
  normalization_id: log_cp10k
```

**metric_configs.yaml**: Contains metric definitions with mathematical descriptions, references, and optimization directions.

**dataset_uns.yaml**: Contains dataset metadata including descriptions, organisms, and source URLs.

### biomlbench Target Structure

Each task should follow this structure:
```
biomlbench/tasks/openproblems/{task_name}/
├── config.yaml
├── leaderboard.csv
├── prepare.py
├── grade.py
├── description.md
└── checksums.yaml (optional)
```

## Implementation Plan

The implementation will be done in the following phases:

1. **Phase 1**: Create OpenProblems data source
2. **Phase 2**: Build data extraction and parsing utilities  
3. **Phase 3**: Implement data format conversion functions
4. **Phase 4**: Create task generation pipeline
5. **Phase 5**: Build ingestion script
6. **Phase 6**: Add validation and testing
7. **Phase 7**: Integrate with biomlbench framework

Each phase will be implemented iteratively with testing at each step.

## Phase 1: Creating the OpenProblems Data Source

### Step 1.1: Create the Base Data Source File

Create `biomlbench/data_sources/openproblems.py`:

```python
"""OpenProblems data source for biomlbench integration."""

import boto3
import yaml
import pandas as pd
from typing import Dict, List, Optional, Tuple
from datetime import datetime
import logging
from .base import BaseDataSource

logger = logging.getLogger(__name__)

class OpenProblemsDataSource(BaseDataSource):
    """Data source for OpenProblems benchmarks stored in AWS S3."""
    
    def __init__(self):
        """Initialize OpenProblems data source."""
        self.bucket_name = "openproblems-data"
        self.base_prefix = "resources"
        self.s3_client = boto3.client('s3', region_name='us-east-1')
        self._task_cache = {}
        
    def discover_tasks(self) -> List[str]:
        """
        Discover all available tasks in OpenProblems S3 bucket.
        
        Returns:
            List of task names (e.g., ['batch_integration', 'label_projection'])
        """
        try:
            response = self.s3_client.list_objects_v2(
                Bucket=self.bucket_name,
                Prefix=f"{self.base_prefix}/",
                Delimiter="/"
            )
            
            tasks = []
            for prefix in response.get('CommonPrefixes', []):
                task_name = prefix['Prefix'].split('/')[-2]
                tasks.append(task_name)
                
            logger.info(f"Discovered {len(tasks)} OpenProblems tasks: {tasks}")
            return sorted(tasks)
            
        except Exception as e:
            logger.error(f"Failed to discover tasks: {e}")
            return []
```

### Step 1.2: Add Run Discovery Methods

Add these methods to the `OpenProblemsDataSource` class:

```python
    def get_latest_run(self, task_name: str) -> Optional[str]:
        """
        Get the latest benchmark run for a specific task.
        
        Args:
            task_name: Name of the task (e.g., 'batch_integration')
            
        Returns:
            S3 prefix for the latest run, or None if not found
        """
        try:
            results_prefix = f"{self.base_prefix}/{task_name}/results/"
            
            response = self.s3_client.list_objects_v2(
                Bucket=self.bucket_name,
                Prefix=results_prefix,
                Delimiter="/"
            )
            
            run_prefixes = []
            for prefix in response.get('CommonPrefixes', []):
                run_dir = prefix['Prefix'].split('/')[-2]
                if run_dir.startswith('run_'):
                    run_prefixes.append(prefix['Prefix'])
            
            if not run_prefixes:
                logger.warning(f"No runs found for task {task_name}")
                return None
                
            # Sort by run timestamp (newest first)
            latest_run = sorted(run_prefixes)[-1]
            logger.info(f"Latest run for {task_name}: {latest_run}")
            return latest_run
            
        except Exception as e:
            logger.error(f"Failed to get latest run for {task_name}: {e}")
            return None
```

### Step 1.3: Add File Download Methods

Add these methods to download and parse S3 files:

```python
    def download_file(self, s3_key: str) -> str:
        """
        Download a file from S3 and return its content as string.
        
        Args:
            s3_key: S3 key for the file
            
        Returns:
            File content as string
        """
        try:
            response = self.s3_client.get_object(
                Bucket=self.bucket_name,
                Key=s3_key
            )
            return response['Body'].read().decode('utf-8')
            
        except Exception as e:
            logger.error(f"Failed to download {s3_key}: {e}")
            raise
    
    def get_task_data(self, task_name: str, run_prefix: Optional[str] = None) -> Dict:
        """
        Download and parse all data files for a task.
        
        Args:
            task_name: Name of the task
            run_prefix: Specific run prefix, or None for latest
            
        Returns:
            Dictionary containing parsed YAML data
        """
        if run_prefix is None:
            run_prefix = self.get_latest_run(task_name)
            if run_prefix is None:
                raise ValueError(f"No runs found for task {task_name}")
        
        # Files to download
        files_to_download = {
            'scores': 'score_uns.yaml',
            'metrics': 'metric_configs.yaml', 
            'methods': 'method_configs.yaml',
            'datasets': 'dataset_uns.yaml',
            'task_info': 'task_info.yaml',
            'state': 'state.yaml'
        }
        
        task_data = {}
        
        for data_type, filename in files_to_download.items():
            s3_key = f"{run_prefix}{filename}"
            
            try:
                content = self.download_file(s3_key)
                task_data[data_type] = yaml.safe_load(content)
                logger.info(f"Downloaded and parsed {filename}")
                
            except Exception as e:
                logger.warning(f"Failed to download {filename}: {e}")
                task_data[data_type] = None
        
        # Add metadata
        task_data['_metadata'] = {
            'task_name': task_name,
            'run_prefix': run_prefix,
            'download_timestamp': datetime.now().isoformat()
        }
        
        return task_data
```

### Step 1.4: Update Data Source Factory

Modify `biomlbench/data_sources/factory.py`:

```python
from .openproblems import OpenProblemsDataSource

def get_data_source(source_name: str):
    """Get a data source instance by name."""
    if source_name == "polaris":
        return PolarisDataSource()
    elif source_name == "openproblems":
        return OpenProblemsDataSource()
    else:
        raise ValueError(f"Unknown data source: {source_name}")

def list_data_sources() -> List[str]:
    """List all available data sources."""
    return ["polaris", "openproblems"]
```

### Step 1.5: Test Phase 1

Create a simple test script to verify the data source works:

```python
#!/usr/bin/env python3
"""Test script for OpenProblems data source."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from biomlbench.data_sources.openproblems import OpenProblemsDataSource

def test_data_source():
    """Test basic functionality of OpenProblems data source."""
    print("Testing OpenProblems data source...")
    
    # Initialize data source
    source = OpenProblemsDataSource()
    
    # Test task discovery
    print("Discovering tasks...")
    tasks = source.discover_tasks()
    print(f"Found {len(tasks)} tasks: {tasks}")
    
    # Test getting latest run for batch_integration
    if 'batch_integration' in tasks:
        print("\nTesting latest run discovery...")
        latest_run = source.get_latest_run('batch_integration')
        print(f"Latest run: {latest_run}")
        
        # Test downloading task data
        print("\nTesting task data download...")
        try:
            task_data = source.get_task_data('batch_integration')
            print(f"Downloaded data keys: {list(task_data.keys())}")
            
            # Check scores
            scores = task_data.get('scores', [])
            print(f"Found {len(scores)} benchmark results")
            
        except Exception as e:
            print(f"Failed to download task data: {e}")
    
    print("Test complete!")

if __name__ == "__main__":
    test_data_source()
```

## Next Steps

After completing Phase 1, test the data source to ensure it can:
1. Discover available OpenProblems tasks
2. Find the latest benchmark run for each task  
3. Download and parse all YAML files

Once Phase 1 is working, proceed to Phase 2: Building the Ingestion Pipeline. 

## Phase 2: Building the Ingestion Pipeline

### Step 2.1: Create Data Conversion Utilities

Create `biomlbench/data_sources/openproblems_utils.py`:

```python
"""Utility functions for OpenProblems data conversion."""

import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Any
import logging

logger = logging.getLogger(__name__)

def extract_metric_names(metrics: List[Dict]) -> List[str]:
    """Extract metric names from metric configs."""
    metric_names = []
    for metric_config in metrics:
        if 'functionality' in metric_config:
            name = metric_config['functionality'].get('name')
            if name:
                metric_names.append(name)
    return sorted(list(set(metric_names)))

def extract_dataset_names(datasets: List[Dict]) -> List[str]:
    """Extract dataset names from dataset info."""
    return sorted([d.get('dataset_id', '') for d in datasets if d.get('dataset_id')])

def aggregate_metrics(metric_ids: List[str], metric_values: List[float], 
                     metric_configs: List[Dict]) -> float:
    """
    Aggregate multiple metric values into a single score.
    
    Args:
        metric_ids: List of metric names
        metric_values: Corresponding metric values  
        metric_configs: Metric configuration data
        
    Returns:
        Aggregated score
    """
    if not metric_values:
        return 0.0
    
    # Create metric metadata lookup
    metric_info = {}
    for config in metric_configs:
        if 'functionality' in config and 'info' in config:
            name = config['functionality'].get('name')
            info_list = config['info'].get('metrics', [])
            if name and info_list:
                info = info_list[0]  # Take first metric info
                metric_info[name] = {
                    'maximize': info.get('maximize', True),
                    'weight': info.get('weight', 1.0),
                    'min_val': info.get('min', 0.0),
                    'max_val': info.get('max', 1.0)
                }
    
    # Normalize and weight metrics
    weighted_scores = []
    total_weight = 0
    
    for metric_id, value in zip(metric_ids, metric_values):
        if metric_id in metric_info:
            info = metric_info[metric_id]
            
            # Normalize to [0, 1] range
            min_val = info['min_val']
            max_val = info['max_val']
            if max_val > min_val:
                normalized = (value - min_val) / (max_val - min_val)
            else:
                normalized = value
            
            # Invert if metric should be minimized
            if not info['maximize']:
                normalized = 1.0 - normalized
            
            # Apply weight
            weight = info['weight']
            weighted_scores.append(normalized * weight)
            total_weight += weight
        else:
            # Unknown metric, use raw value with weight 1
            weighted_scores.append(value)
            total_weight += 1.0
    
    # Return weighted average
    if total_weight > 0:
        return sum(weighted_scores) / total_weight
    else:
        return sum(metric_values) / len(metric_values)

def handle_multi_dataset_scores(scores: List[Dict]) -> Dict[str, List[Dict]]:
    """
    Organize scores by dataset for multi-dataset tasks.
    
    Args:
        scores: List of score entries from OpenProblems
        
    Returns:
        Dictionary mapping dataset_id to list of scores
    """
    dataset_scores = {}
    
    for score in scores:
        dataset_id = score.get('dataset_id', 'unknown')
        if dataset_id not in dataset_scores:
            dataset_scores[dataset_id] = []
        dataset_scores[dataset_id].append(score)
    
    return dataset_scores

def create_leaderboard_dataframe(scores: List[Dict], metrics: List[Dict]) -> pd.DataFrame:
    """
    Create a pandas DataFrame from OpenProblems scores.
    
    Args:
        scores: List of score entries
        metrics: List of metric configurations
        
    Returns:
        DataFrame with columns: teamName, score, submissionDate, dataset, metrics
    """
    if not scores:
        return pd.DataFrame(columns=['teamName', 'score', 'submissionDate'])
    
    leaderboard_rows = []
    
    for score_entry in scores:
        method_id = score_entry.get('method_id', 'unknown')
        date_created = score_entry.get('date_created', '2024-01-01')
        metric_ids = score_entry.get('metric_ids', [])
        metric_values = score_entry.get('metric_values', [])
        dataset_id = score_entry.get('dataset_id', 'unknown')
        
        # Calculate aggregated score
        if len(metric_ids) == len(metric_values) and metric_values:
            agg_score = aggregate_metrics(metric_ids, metric_values, metrics)
        else:
            agg_score = 0.0
        
        leaderboard_rows.append({
            'teamName': method_id,
            'score': agg_score,
            'submissionDate': date_created,
            'dataset': dataset_id,
            'metrics': ','.join(metric_ids),
            'metric_values': ','.join([f"{v:.6f}" for v in metric_values])
        })
    
    df = pd.DataFrame(leaderboard_rows)
    
    # Group by teamName and dataset, take best score
    if len(df) > 0:
        df_agg = df.groupby(['teamName', 'dataset']).agg({
            'score': 'max',
            'submissionDate': 'first',
            'metrics': 'first',
            'metric_values': 'first'
        }).reset_index()
        
        # For multi-dataset tasks, further aggregate by teamName
        final_agg = df_agg.groupby('teamName').agg({
            'score': 'mean',  # Average across datasets
            'submissionDate': 'first',
            'dataset': lambda x: ';'.join(x.unique()),
            'metrics': 'first',
            'metric_values': 'first'
        }).reset_index()
        
        # Sort by score (descending)
        final_agg = final_agg.sort_values('score', ascending=False)
        return final_agg
    
    return df
```

### Step 2.2: Create Task Conversion Class

Create `biomlbench/data_sources/openproblems_converter.py`:

```python
"""Convert OpenProblems data to biomlbench format."""

import yaml
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime
import logging

from .openproblems_utils import (
    extract_metric_names, extract_dataset_names, 
    create_leaderboard_dataframe, handle_multi_dataset_scores
)

logger = logging.getLogger(__name__)

class OpenProblemsConverter:
    """Convert OpenProblems task data to biomlbench format."""
    
    def __init__(self):
        """Initialize converter."""
        pass
    
    def convert_task(self, task_data: Dict) -> Dict[str, Any]:
        """
        Convert OpenProblems task data to biomlbench format.
        
        Args:
            task_data: Raw OpenProblems data
            
        Returns:
            Dictionary with biomlbench-formatted data
        """
        task_name = task_data['_metadata']['task_name']
        logger.info(f"Converting task: {task_name}")
        
        # Extract components
        task_info = task_data.get('task_info', {})
        scores = task_data.get('scores', [])
        metrics = task_data.get('metrics', [])
        datasets = task_data.get('datasets', [])
        methods = task_data.get('methods', [])
        
        # Check if this is a multi-dataset task
        dataset_scores = handle_multi_dataset_scores(scores)
        is_multi_dataset = len(dataset_scores) > 1
        
        if is_multi_dataset:
            logger.info(f"Multi-dataset task with {len(dataset_scores)} datasets")
            return self._convert_multi_dataset_task(task_data, dataset_scores)
        else:
            logger.info("Single dataset task")
            return self._convert_single_task(task_data)
    
    def _convert_single_task(self, task_data: Dict) -> Dict[str, Any]:
        """Convert a single dataset task."""
        task_name = task_data['_metadata']['task_name']
        task_info = task_data.get('task_info', {})
        
        # Generate config.yaml
        config = self._generate_config(task_data)
        
        # Generate leaderboard
        leaderboard_csv = self._generate_leaderboard(task_data)
        
        # Generate scripts
        prepare_script = self._generate_prepare_script(task_data)
        grade_script = self._generate_grade_script(task_data)
        
        # Generate description
        description = self._generate_description(task_data)
        
        return {
            'config.yaml': config,
            'leaderboard.csv': leaderboard_csv,
            'prepare.py': prepare_script,
            'grade.py': grade_script,
            'description.md': description,
            '_metadata': task_data['_metadata']
        }
    
    def _convert_multi_dataset_task(self, task_data: Dict, 
                                  dataset_scores: Dict[str, List[Dict]]) -> Dict[str, Dict]:
        """Convert a multi-dataset task into separate tasks."""
        tasks = {}
        base_task_name = task_data['_metadata']['task_name']
        
        for dataset_id, dataset_score_list in dataset_scores.items():
            # Create clean dataset name
            dataset_name = dataset_id.split('/')[-1].replace('_', '-')
            task_name = f"{base_task_name}-{dataset_name}"
            
            # Create modified task data for this dataset
            modified_task_data = task_data.copy()
            modified_task_data['scores'] = dataset_score_list
            modified_task_data['_metadata'] = task_data['_metadata'].copy()
            modified_task_data['_metadata']['task_name'] = task_name
            modified_task_data['_metadata']['dataset_filter'] = dataset_id
            
            # Convert this dataset-specific task
            tasks[task_name] = self._convert_single_task(modified_task_data)
        
        return tasks
    
    def _generate_config(self, task_data: Dict) -> Dict[str, Any]:
        """Generate config.yaml for a task."""
        task_name = task_data['_metadata']['task_name']
        task_info = task_data.get('task_info', {})
        metrics = task_data.get('metrics', [])
        datasets = task_data.get('datasets', [])
        
        config = {
            'name': f"openproblems-{task_name}",
            'description': task_info.get('summary', f"OpenProblems {task_name} task"),
            'type': 'benchmark',
            'source': 'openproblems',
            'source_url': f"https://openproblems.bio/results/{task_name.replace('-', '_')}",
            'reference': self._extract_references(task_info),
            'data_source': 'openproblems',
            'metrics': extract_metric_names(metrics),
            'datasets': extract_dataset_names(datasets),
            'metadata': {
                'run_prefix': task_data['_metadata']['run_prefix'],
                'download_timestamp': task_data['_metadata']['download_timestamp'],
                'openproblems_task_info': task_info,
                'dataset_filter': task_data['_metadata'].get('dataset_filter')
            }
        }
        
        return config
    
    def _generate_leaderboard(self, task_data: Dict) -> str:
        """Generate leaderboard CSV."""
        scores = task_data.get('scores', [])
        metrics = task_data.get('metrics', [])
        
        if not scores:
            logger.warning("No scores found, creating empty leaderboard")
            return "teamName,score,submissionDate\n"
        
        # Create DataFrame
        df = create_leaderboard_dataframe(scores, metrics)
        
        # Convert to CSV, keeping only essential columns for biomlbench
        if len(df) > 0:
            output_df = df[['teamName', 'score', 'submissionDate']].copy()
            return output_df.to_csv(index=False)
        else:
            return "teamName,score,submissionDate\n"
    
    def _generate_prepare_script(self, task_data: Dict) -> str:
        """Generate prepare.py script."""
        task_name = task_data['_metadata']['task_name']
        datasets = task_data.get('datasets', [])
        
        script = f'''#!/usr/bin/env python3
"""
Data preparation script for OpenProblems {task_name} task.

This script loads and prepares datasets from OpenProblems for evaluation.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

def prepare_data():
    """
    Prepare data for {task_name} task.
    
    Returns:
        Prepared dataset ready for evaluation
    """
    logger.info("Preparing OpenProblems {task_name} dataset")
    
    # Dataset information from OpenProblems
    datasets = {datasets}
    
    logger.info(f"Found {{len(datasets)}} datasets")
    
    # Note: This is a placeholder implementation
    # OpenProblems datasets are stored in AnnData format on S3
    # Full implementation would require downloading and processing
    # the actual benchmark datasets from OpenProblems
    
    return {{
        'datasets': datasets,
        'task_type': '{task_name}',
        'source': 'openproblems',
        'description': 'Placeholder data preparation for OpenProblems task'
    }}

if __name__ == "__main__":
    data = prepare_data()
    print(f"Prepared data for {{data['task_type']}} task")
    print(f"Available datasets: {{len(data['datasets'])}}")
'''
        return script
    
    def _generate_grade_script(self, task_data: Dict) -> str:
        """Generate grade.py script."""
        task_name = task_data['_metadata']['task_name']
        metrics = task_data.get('metrics', [])
        
        # Extract metric information
        metric_info = []
        for metric_config in metrics:
            if 'functionality' in metric_config and 'info' in metric_config:
                func = metric_config['functionality']
                info_list = metric_config['info'].get('metrics', [])
                
                metric_name = func.get('name', 'unknown')
                if info_list:
                    metric_desc = info_list[0]
                    metric_info.append({
                        'name': metric_name,
                        'description': metric_desc.get('description', ''),
                        'maximize': metric_desc.get('maximize', True),
                        'min_val': metric_desc.get('min', 0.0),
                        'max_val': metric_desc.get('max', 1.0)
                    })
        
        script = f'''#!/usr/bin/env python3
"""
Evaluation script for OpenProblems {task_name} task.

This script evaluates method outputs using OpenProblems metrics.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

# Metric definitions from OpenProblems
METRICS = {metric_info}

def grade_submission(predictions, ground_truth):
    """
    Grade a submission for the {task_name} task.
    
    Args:
        predictions: Method predictions
        ground_truth: Ground truth data
        
    Returns:
        Dictionary with metric scores
    """
    logger.info("Grading submission for {task_name}")
    
    # Note: This is a placeholder implementation
    # Full implementation would require implementing the actual
    # OpenProblems metrics and evaluation procedures
    
    scores = {{}}
    
    for metric in METRICS:
        # Calculate score based on metric definition
        min_val = metric.get('min_val', 0.0)
        max_val = metric.get('max_val', 1.0)
        
        # Placeholder: random score within bounds
        scores[metric['name']] = np.random.uniform(min_val, max_val)
        logger.info(f"{{metric['name']}}: {{scores[metric['name']]:.4f}}")
    
    # Calculate overall score (weighted average of all metrics)
    if scores:
        overall_score = np.mean(list(scores.values()))
    else:
        overall_score = 0.0
        
    scores['overall'] = overall_score
    
    return scores

def calculate_medal(score):
    """Calculate medal based on score percentile."""
    if score >= 0.8:
        return "gold"
    elif score >= 0.6:
        return "silver"
    elif score >= 0.4:
        return "bronze"
    else:
        return "none"

if __name__ == "__main__":
    # Example usage
    dummy_predictions = np.random.rand(100, 10)
    dummy_ground_truth = np.random.rand(100, 10)
    
    scores = grade_submission(dummy_predictions, dummy_ground_truth)
    medal = calculate_medal(scores['overall'])
    
    print(f"Overall score: {{scores['overall']:.4f}}")
    print(f"Medal: {{medal}}")
    print(f"Individual metrics: {{{{k: v for k, v in scores.items() if k != 'overall'}}}}")
'''
        return script
    
    def _generate_description(self, task_data: Dict) -> str:
        """Generate description.md file."""
        task_name = task_data['_metadata']['task_name']
        task_info = task_data.get('task_info', {})
        datasets = task_data.get('datasets', [])
        metrics = task_data.get('metrics', [])
        
        # Extract key information
        summary = task_info.get('summary', f'OpenProblems {task_name} task')
        motivation = task_info.get('motivation', 'No motivation provided')
        description = task_info.get('description', 'No description provided')
        
        md_content = f'''# OpenProblems {task_name.replace('_', ' ').replace('-', ' ').title()} Task

## Summary

{summary}

## Motivation

{motivation}

## Description

{description}

## Datasets

This task uses {len(datasets)} datasets from OpenProblems:

'''
        
        # Add dataset information
        for i, dataset in enumerate(datasets[:5]):  # Limit to first 5 for brevity
            dataset_id = dataset.get('dataset_id', 'Unknown')
            dataset_name = dataset.get('dataset_name', 'Unknown')
            dataset_summary = dataset.get('dataset_summary', 'No summary')
            
            md_content += f"{i+1}. **{dataset_name}** (`{dataset_id}`): {dataset_summary}\n"
        
        if len(datasets) > 5:
            md_content += f"\n... and {len(datasets) - 5} more datasets.\n"
        
        # Add metrics information
        md_content += f"\n## Metrics\n\nThis task uses {len(metrics)} evaluation metrics:\n\n"
        
        for i, metric_config in enumerate(metrics[:10]):  # Limit to first 10
            if 'functionality' in metric_config and 'info' in metric_config:
                func = metric_config['functionality']
                info_list = metric_config['info'].get('metrics', [])
                
                metric_name = func.get('name', 'unknown')
                if info_list:
                    metric_desc = info_list[0].get('summary', 'No description')
                    md_content += f"{i+1}. **{metric_name}**: {metric_desc}\n"
        
        # Filter info for dataset-specific tasks
        dataset_filter = task_data['_metadata'].get('dataset_filter')
        if dataset_filter:
            md_content += f"\n**Note**: This task is filtered to dataset `{dataset_filter}`.\n"
        
        md_content += f'''

## Source

This task is automatically generated from OpenProblems benchmarks.

- Original source: https://openproblems.bio/results/{task_name.split('-')[0]}
- Download timestamp: {task_data['_metadata']['download_timestamp']}
- Run: {task_data['_metadata']['run_prefix']}

For more information about OpenProblems, visit https://openproblems.bio/
'''
        
        return md_content
    
    def _extract_references(self, task_info: Dict) -> List[str]:
        """Extract reference information from task info."""
        authors = task_info.get('authors', [])
        if isinstance(authors, list):
            return [author.get('name', 'Unknown') for author in authors if isinstance(author, dict)]
        return []
```

### Step 2.3: Test Phase 2

Create a test script to verify the conversion works:

```python
#!/usr/bin/env python3
"""Test script for OpenProblems data conversion."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from biomlbench.data_sources.openproblems import OpenProblemsDataSource
from biomlbench.data_sources.openproblems_converter import OpenProblemsConverter

def test_conversion():
    """Test OpenProblems data conversion."""
    print("Testing OpenProblems data conversion...")
    
    # Initialize components
    source = OpenProblemsDataSource()
    converter = OpenProblemsConverter()
    
    # Test with batch_integration task
    try:
        print("Downloading batch_integration data...")
        task_data = source.get_task_data('batch_integration')
        
        print("Converting to biomlbench format...")
        converted = converter.convert_task(task_data)
        
        # Check if it's a multi-dataset task
        if isinstance(converted, dict) and 'config.yaml' in converted:
            print("Single task conversion:")
            print(f"  Config keys: {list(converted.keys())}")
            print(f"  Task name: {converted['config.yaml']['name']}")
            print(f"  Metrics: {converted['config.yaml']['metrics']}")
        else:
            print(f"Multi-dataset conversion: {len(converted)} tasks")
            for task_name, task_data in converted.items():
                print(f"  - {task_name}: {list(task_data.keys())}")
        
        print("Conversion test passed!")
        
    except Exception as e:
        print(f"Conversion test failed: {e}")

if __name__ == "__main__":
    test_conversion()
```

## Next Steps for Phase 2

After implementing Phase 2:

1. Test the conversion utilities with different OpenProblems tasks
2. Verify that multi-dataset tasks are handled correctly  
3. Check that leaderboard format matches biomlbench expectations
4. Ensure generated Python scripts have valid syntax

Once Phase 2 is working, proceed to Phase 3: Creating the Main Ingestion Script. 

## Phase 3: Creating the Main Ingestion Script

### Step 3.1: Create the Main Ingestion Script

Create `scripts/ingest_openproblems.py`:

```python
#!/usr/bin/env python3
"""
Ingest OpenProblems benchmarks into biomlbench format.

This script downloads benchmark data from OpenProblems S3 bucket,
converts it to biomlbench format, and creates task directories.
"""

import argparse
import logging
import os
import sys
import yaml
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime

# Add biomlbench to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from biomlbench.data_sources.openproblems import OpenProblemsDataSource
from biomlbench.data_sources.openproblems_converter import OpenProblemsConverter

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class OpenProblemsIngester:
    """Main class for ingesting OpenProblems data."""
    
    def __init__(self, output_dir: str = "biomlbench/tasks/openproblems"):
        """
        Initialize the ingester.
        
        Args:
            output_dir: Directory to create task folders in
        """
        self.output_dir = Path(output_dir)
        self.data_source = OpenProblemsDataSource()
        self.converter = OpenProblemsConverter()
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def process_task(self, task_name: str) -> Dict[str, bool]:
        """
        Process a single OpenProblems task.
        
        Args:
            task_name: Name of the task to process
            
        Returns:
            Dictionary mapping task names to success status
        """
        logger.info(f"Processing task: {task_name}")
        
        try:
            # Download task data
            task_data = self.data_source.get_task_data(task_name)
            
            # Convert to biomlbench format
            converted = self.converter.convert_task(task_data)
            
            # Handle both single and multi-dataset tasks
            results = {}
            
            if isinstance(converted, dict) and 'config.yaml' in converted:
                # Single task
                task_dir = self.output_dir / task_name
                self.write_task_files(task_dir, converted)
                results[task_name] = True
                logger.info(f"✓ Successfully processed {task_name}")
                
            else:
                # Multi-dataset task - multiple sub-tasks
                for sub_task_name, sub_task_data in converted.items():
                    task_dir = self.output_dir / sub_task_name
                    self.write_task_files(task_dir, sub_task_data)
                    results[sub_task_name] = True
                    logger.info(f"✓ Successfully processed {sub_task_name}")
            
            return results
            
        except Exception as e:
            logger.error(f"✗ Failed to process {task_name}: {e}")
            return {task_name: False}
    
    def write_task_files(self, task_dir: Path, task_data: Dict):
        """
        Write all task files to the specified directory.
        
        Args:
            task_dir: Directory to write files to
            task_data: Converted task data
        """
        task_dir.mkdir(parents=True, exist_ok=True)
        
        # Write YAML config
        config_file = task_dir / 'config.yaml'
        with open(config_file, 'w') as f:
            yaml.dump(task_data['config.yaml'], f, default_flow_style=False)
        
        # Write CSV leaderboard
        leaderboard_file = task_dir / 'leaderboard.csv'
        with open(leaderboard_file, 'w') as f:
            f.write(task_data['leaderboard.csv'])
        
        # Write Python scripts
        for script_name in ['prepare.py', 'grade.py']:
            script_file = task_dir / script_name
            with open(script_file, 'w') as f:
                f.write(task_data[script_name])
            # Make executable
            script_file.chmod(0o755)
        
        # Write description
        desc_file = task_dir / 'description.md'
        with open(desc_file, 'w') as f:
            f.write(task_data['description.md'])
        
        logger.info(f"Wrote task files to {task_dir}")
    
    def run_full_ingestion(self, tasks: Optional[List[str]] = None) -> Dict[str, bool]:
        """
        Run full ingestion for all or specified tasks.
        
        Args:
            tasks: List of specific tasks to process, or None for all
            
        Returns:
            Dictionary mapping task names to success status
        """
        if tasks is None:
            tasks = self.data_source.discover_tasks()
        
        logger.info(f"Starting ingestion for {len(tasks)} tasks")
        
        all_results = {}
        successful_tasks = 0
        total_tasks = 0
        
        for task_name in tasks:
            task_results = self.process_task(task_name)
            all_results.update(task_results)
            
            # Count successes
            task_successes = sum(task_results.values())
            successful_tasks += task_successes
            total_tasks += len(task_results)
        
        logger.info(f"Ingestion complete: {successful_tasks}/{total_tasks} tasks successful")
        
        # Write summary
        self.write_ingestion_summary(all_results, successful_tasks, total_tasks)
        
        return all_results
    
    def write_ingestion_summary(self, results: Dict[str, bool], 
                              successful: int, total: int):
        """Write ingestion summary to file."""
        summary_file = self.output_dir / 'ingestion_summary.yaml'
        summary = {
            'timestamp': datetime.now().isoformat(),
            'total_tasks': total,
            'successful_tasks': successful,
            'failed_tasks': total - successful,
            'success_rate': f"{successful/total*100:.1f}%" if total > 0 else "0%",
            'results': results,
            'successful_task_names': [name for name, success in results.items() if success],
            'failed_task_names': [name for name, success in results.items() if not success]
        }
        
        with open(summary_file, 'w') as f:
            yaml.dump(summary, f, default_flow_style=False)
        
        logger.info(f"Wrote ingestion summary to {summary_file}")

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description='Ingest OpenProblems benchmarks')
    parser.add_argument('--output-dir', default='biomlbench/tasks/openproblems',
                       help='Output directory for tasks')
    parser.add_argument('--tasks', multiple=True,
                       help='Specific tasks to ingest (default: all)')
    parser.add_argument('--update', is_flag=True,
                       help='Update existing tasks')
    parser.add_argument('--validate', is_flag=True,
                       help='Validate after ingestion')
    parser.add_argument('--verbose', '-v', action='store_true')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    ingester = OpenProblemsIngester(args.output_dir)
    
    if args.update:
        # Update existing tasks
        logger.info("Updating existing tasks...")
        results = ingester.run_full_ingestion()
    else:
        # Ingest new tasks
        logger.info("Ingesting new tasks...")
        results = ingester.run_full_ingestion(args.tasks)
    
    # Print summary
    successful = sum(results.values())
    total = len(results)
    
    print(f"\nIngestion Summary:")
    print(f"  Total tasks: {total}")
    print(f"  Successful: {successful}")
    print(f"  Failed: {total - successful}")
    print(f"  Success rate: {successful/total*100:.1f}%" if total > 0 else "  Success rate: 0%")
    
    if total - successful > 0:
        print(f"\nFailed tasks:")
        for task, success in results.items():
            if not success:
                print(f"  - {task}")
    
    # Validate after ingestion if requested
    if args.validate:
        logger.info("Running validation...")
        validate_results = subprocess.run([sys.executable, 'scripts/validate_openproblems_ingestion.py', args.output_dir], capture_output=False)
        return validate_results.returncode
    
    return 0 if successful == total else 1

if __name__ == "__main__":
    sys.exit(main())
```

### Step 3.2: Create Helper Scripts

Create `scripts/update_openproblems.py` for updating existing tasks:

```python
#!/usr/bin/env python3
"""
Update existing OpenProblems tasks with latest data.

This script checks for newer benchmark runs and updates tasks if needed.
"""

import argparse
import logging
import sys
import yaml
from pathlib import Path
from datetime import datetime

# Add biomlbench to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.ingest_openproblems import OpenProblemsIngester

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def check_for_updates(task_dir: Path, data_source) -> bool:
    """
    Check if a task needs updating.
    
    Args:
        task_dir: Path to existing task directory
        data_source: OpenProblems data source
        
    Returns:
        True if task needs updating
    """
    config_file = task_dir / 'config.yaml'
    if not config_file.exists():
        return True  # No config, needs update
    
    try:
        with open(config_file) as f:
            config = yaml.safe_load(f)
        
        # Get stored run prefix
        stored_run = config.get('metadata', {}).get('run_prefix', '')
        
        # Get latest run
        task_name = task_dir.name
        if task_name.startswith('openproblems-'):
            base_task_name = task_name.replace('openproblems-', '').split('-')[0]
        else:
            base_task_name = task_name.split('-')[0]
        
        latest_run = data_source.get_latest_run(base_task_name)
        
        if latest_run and latest_run != stored_run:
            logger.info(f"Task {task_name} has newer data available")
            logger.info(f"  Current: {stored_run}")
            logger.info(f"  Latest:  {latest_run}")
            return True
        
        return False
        
    except Exception as e:
        logger.warning(f"Failed to check update status for {task_name}: {e}")
        return True  # If we can't check, assume update needed

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description='Update OpenProblems tasks')
    parser.add_argument('--task-dir', default='biomlbench/tasks/openproblems',
                       help='Directory containing OpenProblems tasks')
    parser.add_argument('--force', action='store_true',
                       help='Force update all tasks')
    parser.add_argument('--verbose', '-v', action='store_true')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    task_dir = Path(args.task_dir)
    if not task_dir.exists():
        logger.error(f"Task directory does not exist: {task_dir}")
        return 1
    
    # Initialize ingester
    ingester = OpenProblemsIngester(str(task_dir))
    
    # Find existing tasks
    existing_tasks = [d for d in task_dir.iterdir() if d.is_dir()]
    logger.info(f"Found {len(existing_tasks)} existing tasks")
    
    # Check which tasks need updates
    tasks_to_update = []
    
    for task_path in existing_tasks:
        if args.force or check_for_updates(task_path, ingester.data_source):
            tasks_to_update.append(task_path.name)
    
    if not tasks_to_update:
        logger.info("No tasks need updating")
        return 0
    
    logger.info(f"Updating {len(tasks_to_update)} tasks: {tasks_to_update}")
    
    # Extract base task names for updating
    base_tasks = set()
    for task_name in tasks_to_update:
        if task_name.startswith('openproblems-'):
            base_name = task_name.replace('openproblems-', '').split('-')[0]
        else:
            base_name = task_name.split('-')[0]
        base_tasks.add(base_name)
    
    # Run updates
    results = ingester.run_full_ingestion(list(base_tasks))
    
    # Report results
    successful = sum(results.values())
    total = len(results)
    
    print(f"\nUpdate Summary:")
    print(f"  Tasks checked: {len(existing_tasks)}")
    print(f"  Tasks updated: {successful}")
    print(f"  Update failures: {total - successful}")
    
    return 0 if successful == total else 1

if __name__ == "__main__":
    sys.exit(main())
```

## Phase 4: Testing and Validation

### Step 4.1: Create Comprehensive Tests

Create `tests/test_openproblems_integration.py`:

```python
"""Tests for OpenProblems integration."""

import pytest
import tempfile
import yaml
import pandas as pd
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

# Import the classes we're testing
from biomlbench.data_sources.openproblems import OpenProblemsDataSource
from biomlbench.data_sources.openproblems_converter import OpenProblemsConverter
from biomlbench.data_sources.openproblems_utils import (
    extract_metric_names, aggregate_metrics, handle_multi_dataset_scores
)
from scripts.ingest_openproblems import OpenProblemsIngester

class TestOpenProblemsDataSource:
    """Test OpenProblems data source functionality."""
    
    @patch('boto3.client')
    def test_discover_tasks(self, mock_boto):
        """Test task discovery."""
        mock_client = Mock()
        mock_boto.return_value = mock_client
        
        mock_client.list_objects_v2.return_value = {
            'CommonPrefixes': [
                {'Prefix': 'resources/batch_integration/'},
                {'Prefix': 'resources/label_projection/'},
                {'Prefix': 'resources/perturbation_prediction/'}
            ]
        }
        
        source = OpenProblemsDataSource()
        tasks = source.discover_tasks()
        
        assert tasks == ['batch_integration', 'label_projection', 'perturbation_prediction']
        mock_client.list_objects_v2.assert_called_once()
    
    @patch('boto3.client')
    def test_get_latest_run(self, mock_boto):
        """Test getting latest run."""
        mock_client = Mock()
        mock_boto.return_value = mock_client
        
        mock_client.list_objects_v2.return_value = {
            'CommonPrefixes': [
                {'Prefix': 'resources/batch_integration/results/run_2024-06-27_10-30-15/'},
                {'Prefix': 'resources/batch_integration/results/run_2024-06-28_13-20-27/'},
                {'Prefix': 'resources/batch_integration/results/run_2024-06-26_08-15-42/'}
            ]
        }
        
        source = OpenProblemsDataSource()
        latest = source.get_latest_run('batch_integration')
        
        assert latest == 'resources/batch_integration/results/run_2024-06-28_13-20-27/'

class TestUtilityFunctions:
    """Test utility functions."""
    
    def test_extract_metric_names(self):
        """Test extracting metric names from configs."""
        metrics = [
            {'functionality': {'name': 'metric_a'}},
            {'functionality': {'name': 'metric_b'}},
            {'other_data': 'value'}  # Should be ignored
        ]
        
        names = extract_metric_names(metrics)
        assert names == ['metric_a', 'metric_b']
    
    def test_aggregate_metrics(self):
        """Test metric aggregation."""
        metric_ids = ['metric_a', 'metric_b']
        metric_values = [0.8, 0.6]
        metric_configs = [
            {
                'functionality': {'name': 'metric_a'},
                'info': {'metrics': [{'maximize': True, 'min': 0, 'max': 1}]}
            },
            {
                'functionality': {'name': 'metric_b'},
                'info': {'metrics': [{'maximize': True, 'min': 0, 'max': 1}]}
            }
        ]
        
        score = aggregate_metrics(metric_ids, metric_values, metric_configs)
        assert 0.0 <= score <= 1.0
        assert abs(score - 0.7) < 0.1  # Should be approximately average
    
    def test_handle_multi_dataset_scores(self):
        """Test multi-dataset score handling."""
        scores = [
            {'dataset_id': 'dataset_a', 'method_id': 'method_1', 'score': 0.8},
            {'dataset_id': 'dataset_b', 'method_id': 'method_1', 'score': 0.6},
            {'dataset_id': 'dataset_a', 'method_id': 'method_2', 'score': 0.7}
        ]
        
        dataset_scores = handle_multi_dataset_scores(scores)
        
        assert 'dataset_a' in dataset_scores
        assert 'dataset_b' in dataset_scores
        assert len(dataset_scores['dataset_a']) == 2
        assert len(dataset_scores['dataset_b']) == 1

class TestDataConversion:
    """Test data conversion functionality."""
    
    def test_single_task_conversion(self):
        """Test converting a single-dataset task."""
        converter = OpenProblemsConverter()
        
        # Sample data
        task_data = {
            'scores': [
                {
                    'dataset_id': 'test_dataset',
                    'method_id': 'test_method',
                    'metric_ids': ['metric_a'],
                    'metric_values': [0.8],
                    'date_created': '2024-01-01'
                }
            ],
            'task_info': {'summary': 'Test task'},
            'metrics': [
                {
                    'functionality': {'name': 'metric_a'},
                    'info': {'metrics': [{'maximize': True, 'min': 0, 'max': 1}]}
                }
            ],
            'datasets': [{'dataset_id': 'test_dataset'}],
            'methods': [],
            '_metadata': {
                'task_name': 'test_task',
                'run_prefix': 'test/',
                'download_timestamp': '2024-01-01T00:00:00'
            }
        }
        
        result = converter.convert_task(task_data)
        
        # Should be single task (not multi-dataset)
        assert 'config.yaml' in result
        assert 'leaderboard.csv' in result
        assert 'prepare.py' in result
        assert 'grade.py' in result
        assert 'description.md' in result
        
        # Check config content
        config = result['config.yaml']
        assert config['name'] == 'openproblems-test_task'
        assert config['source'] == 'openproblems'
    
    def test_multi_dataset_conversion(self):
        """Test converting a multi-dataset task."""
        converter = OpenProblemsConverter()
        
        # Sample data with multiple datasets
        task_data = {
            'scores': [
                {
                    'dataset_id': 'dataset_a',
                    'method_id': 'method_1',
                    'metric_ids': ['metric_a'],
                    'metric_values': [0.8],
                    'date_created': '2024-01-01'
                },
                {
                    'dataset_id': 'dataset_b',
                    'method_id': 'method_1',
                    'metric_ids': ['metric_a'],
                    'metric_values': [0.6],
                    'date_created': '2024-01-01'
                }
            ],
            'task_info': {'summary': 'Test multi-dataset task'},
            'metrics': [],
            'datasets': [],
            'methods': [],
            '_metadata': {
                'task_name': 'test_task',
                'run_prefix': 'test/',
                'download_timestamp': '2024-01-01T00:00:00'
            }
        }
        
        result = converter.convert_task(task_data)
        
        # Should be multi-task (dictionary of tasks)
        assert isinstance(result, dict)
        assert 'test_task-dataset-a' in result or 'test_task-dataset_a' in result
        assert len(result) == 2  # Two datasets

class TestIngestionPipeline:
    """Test full ingestion pipeline."""
    
    def test_file_writing(self):
        """Test that files are written correctly."""
        with tempfile.TemporaryDirectory() as temp_dir:
            ingester = OpenProblemsIngester(temp_dir)
            
            # Sample task data
            task_data = {
                'config.yaml': {'name': 'test-task', 'description': 'Test'},
                'leaderboard.csv': 'teamName,score,submissionDate\ntest,0.8,2024-01-01\n',
                'prepare.py': '#!/usr/bin/env python3\nprint("prepare")\n',
                'grade.py': '#!/usr/bin/env python3\nprint("grade")\n',
                'description.md': '# Test Task\n\nTest description\n'
            }
            
            task_dir = Path(temp_dir) / 'test_task'
            ingester.write_task_files(task_dir, task_data)
            
            # Check that files were created
            assert (task_dir / 'config.yaml').exists()
            assert (task_dir / 'leaderboard.csv').exists()
            assert (task_dir / 'prepare.py').exists()
            assert (task_dir / 'grade.py').exists()
            assert (task_dir / 'description.md').exists()
            
            # Check that scripts are executable
            assert (task_dir / 'prepare.py').stat().st_mode & 0o755
            assert (task_dir / 'grade.py').stat().st_mode & 0o755
            
            # Check content
            with open(task_dir / 'config.yaml') as f:
                config = yaml.safe_load(f)
                assert config['name'] == 'test-task'

def test_script_syntax():
    """Test that generated scripts have valid Python syntax."""
    converter = OpenProblemsConverter()
    
    task_data = {
        'task_info': {'summary': 'Test task'},
        'metrics': [],
        'datasets': [],
        '_metadata': {
            'task_name': 'test_task',
            'run_prefix': 'test/',
            'download_timestamp': '2024-01-01T00:00:00'
        }
    }
    
    prepare_script = converter._generate_prepare_script(task_data)
    grade_script = converter._generate_grade_script(task_data)
    
    # Check syntax
    compile(prepare_script, '<prepare.py>', 'exec')
    compile(grade_script, '<grade.py>', 'exec')

if __name__ == "__main__":
    pytest.main([__file__])
```

### Step 4.2: Create Validation Script

Create `scripts/validate_openproblems_ingestion.py`:

```python
#!/usr/bin/env python3
"""
Validate OpenProblems ingestion results.

This script checks that ingested tasks are properly formatted
and compatible with biomlbench.
"""

import argparse
import logging
import yaml
import pandas as pd
from pathlib import Path
import sys
import importlib.util

# Add biomlbench to path
sys.path.insert(0, str(Path(__file__).parent.parent))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def validate_task_structure(task_dir: Path) -> bool:
    """Validate that a task directory has the correct structure."""
    required_files = ['config.yaml', 'leaderboard.csv', 'prepare.py', 'grade.py']
    
    logger.info(f"Validating task structure: {task_dir.name}")
    
    valid = True
    
    # Check required files
    for filename in required_files:
        filepath = task_dir / filename
        if not filepath.exists():
            logger.error(f"Missing required file: {filename}")
            valid = False
    
    return valid

def validate_config_format(task_dir: Path) -> bool:
    """Validate config.yaml format."""
    config_file = task_dir / 'config.yaml'
    if not config_file.exists():
        return False
    
    try:
        with open(config_file) as f:
            config = yaml.safe_load(f)
        
        required_keys = ['name', 'description', 'source']
        for key in required_keys:
            if key not in config:
                logger.error(f"Missing required config key: {key}")
                return False
        
        # Check OpenProblems-specific requirements
        if config.get('source') != 'openproblems':
            logger.warning(f"Expected source='openproblems', got {config.get('source')}")
        
        return True
        
    except Exception as e:
        logger.error(f"Invalid config.yaml: {e}")
        return False

def validate_leaderboard_format(task_dir: Path) -> bool:
    """Validate leaderboard.csv format."""
    leaderboard_file = task_dir / 'leaderboard.csv'
    if not leaderboard_file.exists():
        return False
    
    try:
        df = pd.read_csv(leaderboard_file)
        required_columns = ['teamName', 'score', 'submissionDate']
        
        for col in required_columns:
            if col not in df.columns:
                logger.error(f"Missing required leaderboard column: {col}")
                return False
        
        # Check data types and ranges
        if len(df) > 0:
            if not pd.api.types.is_numeric_dtype(df['score']):
                logger.error("Score column must be numeric")
                return False
            
            # Check score range
            if df['score'].min() < 0 or df['score'].max() > 1:
                logger.warning(f"Scores outside [0,1] range: {df['score'].min():.3f} to {df['score'].max():.3f}")
        
        logger.info(f"✓ Leaderboard has {len(df)} entries")
        return True
        
    except Exception as e:
        logger.error(f"Invalid leaderboard.csv: {e}")
        return False

def validate_script_syntax(task_dir: Path) -> bool:
    """Validate Python script syntax."""
    valid = True
    
    for script_name in ['prepare.py', 'grade.py']:
        script_file = task_dir / script_name
        if script_file.exists():
            try:
                with open(script_file) as f:
                    compile(f.read(), script_file, 'exec')
                logger.debug(f"✓ {script_name} has valid syntax")
            except SyntaxError as e:
                logger.error(f"Syntax error in {script_name}: {e}")
                valid = False
        else:
            logger.error(f"Missing script: {script_name}")
            valid = False
    
    return valid

def validate_task_complete(task_dir: Path) -> bool:
    """Run all validations for a task."""
    validations = [
        validate_task_structure,
        validate_config_format,
        validate_leaderboard_format,
        validate_script_syntax
    ]
    
    all_valid = True
    for validation_func in validations:
        if not validation_func(task_dir):
            all_valid = False
    
    return all_valid

def main():
    """Main validation script."""
    parser = argparse.ArgumentParser(description='Validate OpenProblems ingestion')
    parser.add_argument('task_dir', nargs='?', 
                       default='biomlbench/tasks/openproblems',
                       help='Directory containing ingested tasks')
    parser.add_argument('--task', help='Validate specific task only')
    parser.add_argument('--verbose', '-v', action='store_true')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    task_dir = Path(args.task_dir)
    
    if not task_dir.exists():
        logger.error(f"Task directory does not exist: {task_dir}")
        return 1
    
    # Find tasks to validate
    if args.task:
        tasks_to_validate = [task_dir / args.task]
    else:
        tasks_to_validate = [d for d in task_dir.iterdir() if d.is_dir()]
    
    logger.info(f"Validating {len(tasks_to_validate)} tasks")
    
    # Validate each task
    results = {}
    for task_path in tasks_to_validate:
        task_name = task_path.name
        logger.info(f"\nValidating task: {task_name}")
        
        is_valid = validate_task_complete(task_path)
        results[task_name] = is_valid
        
        if is_valid:
            logger.info(f"✓ {task_name} validation passed")
        else:
            logger.error(f"✗ {task_name} validation failed")
    
    # Summary
    total_tasks = len(results)
    valid_tasks = sum(results.values())
    
    print(f"\nValidation Summary:")
    print(f"  Total tasks: {total_tasks}")
    print(f"  Valid tasks: {valid_tasks}")
    print(f"  Invalid tasks: {total_tasks - valid_tasks}")
    
    if total_tasks - valid_tasks > 0:
        print(f"\nInvalid tasks:")
        for task_name, is_valid in results.items():
            if not is_valid:
                print(f"  - {task_name}")
    
    return 0 if valid_tasks == total_tasks else 1

if __name__ == "__main__":
    sys.exit(main())
```

## Phase 5: Final Integration

### Step 5.1: Update biomlbench Registry

Update `biomlbench/data_sources/__init__.py`:

```python
"""Data sources for biomlbench."""

from .factory import get_data_source, list_data_sources
from .polaris import PolarisDataSource
from .openproblems import OpenProblemsDataSource

__all__ = [
    'get_data_source',
    'list_data_sources', 
    'PolarisDataSource',
    'OpenProblemsDataSource'
]
```

### Step 5.2: Add CLI Integration

Update `biomlbench/cli.py` to include OpenProblems commands:

```python
# Add to existing CLI commands

@click.command()
@click.option('--output-dir', default='biomlbench/tasks/openproblems',
              help='Output directory for tasks')
@click.option('--tasks', multiple=True,
              help='Specific tasks to ingest (default: all)')
@click.option('--update', is_flag=True,
              help='Update existing tasks')
@click.option('--validate', is_flag=True,
              help='Validate after ingestion')
@click.option('--verbose', '-v', is_flag=True,
              help='Verbose output')
def ingest_openproblems(output_dir, tasks, update, validate, verbose):
    """Ingest OpenProblems benchmarks."""
    import subprocess
    import sys
    
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Prepare command
    cmd = [sys.executable, 'scripts/ingest_openproblems.py']
    cmd.extend(['--output-dir', output_dir])
    
    if tasks:
        cmd.extend(['--tasks'] + list(tasks))
    
    if verbose:
        cmd.append('--verbose')
    
    # Run ingestion
    logger.info("Running OpenProblems ingestion...")
    result = subprocess.run(cmd, capture_output=False)
    
    if result.returncode != 0:
        logger.error("Ingestion failed")
        return result.returncode
    
    # Run validation if requested
    if validate:
        logger.info("Running validation...")
        validate_cmd = [sys.executable, 'scripts/validate_openproblems_ingestion.py', output_dir]
        if verbose:
            validate_cmd.append('--verbose')
        
        validate_result = subprocess.run(validate_cmd, capture_output=False)
        return validate_result.returncode
    
    return 0

# Add to main CLI group
cli.add_command(ingest_openproblems)
```

### Step 5.3: Update Documentation

Create `docs/data_sources/openproblems.md`:

```markdown
# OpenProblems Data Source

The OpenProblems data source provides access to benchmark results from the OpenProblems project, a community-driven platform for benchmarking single-cell analysis methods.

## Overview

OpenProblems stores comprehensive benchmark data including:
- Method performance scores across multiple metrics
- Dataset metadata and descriptions  
- Metric definitions with mathematical formulations
- Method configurations and dependencies

## Usage

### Ingesting OpenProblems Tasks

To ingest all available OpenProblems tasks:

```bash
biomlbench ingest-openproblems
```

To ingest specific tasks:

```bash
biomlbench ingest-openproblems --tasks batch_integration label_projection
```

### Available Tasks

OpenProblems includes tasks such as:
- **batch_integration**: Remove batch effects while preserving biological variation
- **label_projection**: Predict cell type labels for new datasets
- **perturbation_prediction**: Predict effects of genetic perturbations
- **dimensionality_reduction**: Reduce data dimensionality while preserving structure

### Data Structure

Each ingested task includes:
- `config.yaml`: Task configuration and metadata
- `leaderboard.csv`: Method performance rankings
- `prepare.py`: Data preparation script
- `grade.py`: Evaluation script
- `description.md`: Task description and documentation

## Implementation Details

The OpenProblems integration:
1. Discovers available tasks from S3 storage
2. Downloads benchmark results and metadata
3. Converts multi-metric scores to single aggregated scores
4. Handles multi-dataset tasks by creating separate task instances
5. Generates biomlbench-compatible task structures

## Updating Tasks

To update existing tasks with latest benchmark data:

```bash
python scripts/update_openproblems.py
```

## Validation

To validate ingested tasks:

```bash
python scripts/validate_openproblems_ingestion.py biomlbench/tasks/openproblems
```
```

## Final Summary

This comprehensive guide provides everything needed to implement OpenProblems integration:

1. **Phase 1**: OpenProblems data source for S3 access
2. **Phase 2**: Data conversion utilities and task converter
3. **Phase 3**: Main ingestion script and helper tools
4. **Phase 4**: Comprehensive testing and validation
5. **Phase 5**: Integration with biomlbench framework

The implementation handles:
- ✅ Multi-metric aggregation
- ✅ Multi-dataset task splitting
- ✅ Error handling and logging
- ✅ Validation and testing
- ✅ CLI integration
- ✅ Update mechanisms

To complete the implementation:
1. Implement each phase sequentially
2. Test each phase before proceeding
3. Run validation after ingestion
4. Update documentation as needed