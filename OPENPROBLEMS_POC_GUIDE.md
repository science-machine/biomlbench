# OpenProblems Proof of Concept Integration Guide

## Goal
Implement **ONE** OpenProblems task (`batch_integration`) as a results-only benchmark in biomlbench. This is a proof of concept focusing on simplicity and speed.

## Approach
- **Results-only integration**: Extract benchmark results from OpenProblems S3, display them in biomlbench
- **No executable tasks**: Don't try to recreate OpenProblems' complex Viash/AnnData pipeline
- **Single task focus**: Just `batch_integration` to prove the concept works
- **Manual implementation**: No automation, no testing, minimal error handling

## Architecture Overview

```
biomlbench/tasks/openproblems/batch_integration/
├── config.yaml          # Task metadata and info
├── leaderboard.csv       # Results from OpenProblems S3
├── description.md        # Task description
├── prepare.py           # Data preparation script
├── grade.py             # Grading script
└── results_data.yaml    # Raw OpenProblems data for reference
```

## Implementation Steps

### Step 1: Create OpenProblems Data Source

Create `biomlbench/data_sources/openproblems.py`:

```python
"""Simple OpenProblems data source for PoC."""

import boto3
import yaml
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Any
import logging

from biomlbench.data_sources.base import DataSource

logger = logging.getLogger(__name__)

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
        task_name = source_config.get('task_name', 'batch_integration')
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
            try:
                response = self.s3_client.get_object(Bucket=self.bucket_name, Key=s3_key)
                content = response['Body'].read().decode('utf-8')
                task_data[data_type] = yaml.safe_load(content)
                logger.info(f"Downloaded {filename}")
            except Exception as e:
                logger.error(f"Failed to download {filename}: {e}")
                task_data[data_type] = []
        
        return task_data
    
    def create_simple_leaderboard(self, task_data: Dict) -> pd.DataFrame:
        """Create a simple leaderboard preserving individual metrics."""
        scores = task_data.get('scores', [])
        
        # Create one row per method-dataset-metric combination
        rows = []
        for entry in scores:
            method_id = entry.get('method_id', 'unknown')
            dataset_id = entry.get('dataset_id', 'unknown')
            metric_ids = entry.get('metric_ids', [])
            metric_values = entry.get('metric_values', [])
            date_created = entry.get('date_created', '2024-06-28')
            
            # Create row for each metric
            for metric_id, metric_value in zip(metric_ids, metric_values):
                rows.append({
                    'method': method_id,
                    'dataset': dataset_id, 
                    'metric': metric_id,
                    'score': metric_value,
                    'date': date_created
                })
        
        return pd.DataFrame(rows)
```

### Step 2: Register Data Source

Add the OpenProblems data source to the factory in `biomlbench/data_sources/factory.py`:

```python
from .openproblems import OpenProblemsDataSource

@register_data_source("openproblems")
class OpenProblemsDataSource(DataSource):
    # ... (class definition from Step 1)
```

Or if there's no decorator system, manually add it to the factory function.

### Step 3: Create Task Generation Script

Create `scripts/create_openproblems_poc.py`:

```python
#!/usr/bin/env python3
"""Create OpenProblems batch_integration task for PoC."""

import yaml
import pandas as pd
from pathlib import Path
import sys

# Add biomlbench to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from biomlbench.data_sources.openproblems import OpenProblemsDataSource

def create_batch_integration_task():
    """Create the batch_integration task."""
    
    # Create task directory
    task_dir = Path("biomlbench/tasks/openproblems/batch_integration")
    task_dir.mkdir(parents=True, exist_ok=True)
    
    # Download data
    print("Downloading OpenProblems batch_integration data...")
    source = OpenProblemsDataSource()
    task_data = source.download_task_results("batch_integration")
    
    # Create leaderboard
    print("Creating leaderboard...")
    leaderboard_df = source.create_simple_leaderboard(task_data)
    
    # Save raw data for reference
    with open(task_dir / "results_data.yaml", "w") as f:
        yaml.dump(task_data, f, default_flow_style=False)
    
    # Create simplified leaderboard CSV (biomlbench format)
    # Group by method and take first score of first metric for simplicity
    simple_leaderboard = []
    for method in leaderboard_df['method'].unique():
        method_data = leaderboard_df[leaderboard_df['method'] == method]
        if len(method_data) > 0:
            # Just take the first score for PoC
            first_score = method_data.iloc[0]
            simple_leaderboard.append({
                'teamName': method,
                'score': first_score['score'],
                'submissionDate': first_score['date']
            })
    
    simple_df = pd.DataFrame(simple_leaderboard)
    simple_df.to_csv(task_dir / "leaderboard.csv", index=False)
    
    # Create config.yaml
    task_info = task_data.get('task_info', {})
    config = {
        'id': 'openproblems/batch_integration',
        'name': 'OpenProblems Batch Integration',
        'task_type': 'benchmark_results',
        'domain': 'single_cell',
        'difficulty': 'reference',
        'description': task_info.get('summary', 'Batch integration benchmark from OpenProblems'),
        'dataset': {
            'answers': 'openproblems/batch_integration/prepared/private/answers.csv',
            'sample_submission': 'openproblems/batch_integration/prepared/public/sample_submission.csv'
        },
        'data_source': {
            'type': 'openproblems',
            'task_name': 'batch_integration'
        },
        'grader': {
            'name': 'reference-only',
            'grade_fn': 'biomlbench.tasks.openproblems.batch_integration.grade:grade'
        },
        'preparer': 'biomlbench.tasks.openproblems.batch_integration.prepare:prepare',
        'biomedical_metadata': {
            'modality': 'single_cell',
            'data_type': 'benchmark_results',
            'source': 'openproblems'
        }
    }
    
    with open(task_dir / "config.yaml", "w") as f:
        yaml.dump(config, f, default_flow_style=False)
    
    # Create description.md
    task_info = task_data.get('task_info', {})
    description = f"""# OpenProblems Batch Integration

## Summary
{task_info.get('summary', 'Batch integration benchmark from OpenProblems')}

## Description  
{task_info.get('description', 'This task evaluates batch integration methods on their ability to remove batch effects while preserving biological variation.')}

## Motivation
{task_info.get('motivation', 'Single-cell datasets often contain batch effects that must be computationally removed.')}

## Data Source
This is a **results-only** benchmark extracted from OpenProblems. The leaderboard shows method performance from the OpenProblems evaluation.

**Note**: This task is for viewing benchmark results only. It cannot be executed by agents.

## Metrics
The full results include multiple metrics per method-dataset combination:
- asw_batch: Average silhouette width for batches
- asw_label: Average silhouette width for labels  
- graph_connectivity: Graph connectivity preservation
- kbet: k-nearest neighbor batch effect test
- And others...

## Source
- Original: https://openproblems.bio/results/batch_integration
- Data: s3://openproblems-data/resources/batch_integration/
- Run: 2024-06-28_13-20-27

For full details and to run evaluations, visit OpenProblems.bio.
"""
    
    with open(task_dir / "description.md", "w") as f:
        f.write(description)
    
    # Create minimal prepare.py (reference-only)
    prepare_script = '''#!/usr/bin/env python3
"""
Reference-only data preparation for OpenProblems batch_integration.

This task is for viewing benchmark results only.
For actual data and evaluation, use OpenProblems directly.
"""

from pathlib import Path
import pandas as pd
import logging

logger = logging.getLogger(__name__)

def prepare(raw: Path, public: Path, private: Path) -> None:
    """Prepare reference data (no actual data preparation)."""
    logger.info("This is a reference-only task from OpenProblems")
    logger.info("No data preparation needed - results are pre-computed")
    logger.info("Visit https://openproblems.bio for actual evaluation")
    
    # Create directories
    public.mkdir(parents=True, exist_ok=True)
    private.mkdir(parents=True, exist_ok=True)
    
    # Create dummy sample_submission.csv
    dummy_df = pd.DataFrame({'id': [1], 'prediction': [0.5]})
    dummy_df.to_csv(public / "sample_submission.csv", index=False)
    
    # Create dummy answers.csv  
    dummy_df.to_csv(private / "answers.csv", index=False)

if __name__ == "__main__":
    print("OpenProblems Batch Integration - Reference Only")
    print("This task displays pre-computed benchmark results")
'''
    
    with open(task_dir / "prepare.py", "w") as f:
        f.write(prepare_script)
    
    # Create minimal grade.py (reference-only)  
    grade_script = '''#!/usr/bin/env python3
"""
Reference-only grading for OpenProblems batch_integration.

This task is for viewing benchmark results only.
"""

import pandas as pd
import logging
from typing import Any

logger = logging.getLogger(__name__)

def grade(submission: pd.DataFrame, answers: Any) -> float:
    """Reference-only grading (always returns 0)."""
    logger.info("This is a reference-only task from OpenProblems")
    logger.info("No grading available - results are pre-computed")
    logger.info("Visit https://openproblems.bio for actual evaluation")
    
    return 0.0

if __name__ == "__main__":
    print("OpenProblems Batch Integration - Reference Only")
    print("This task displays pre-computed benchmark results")
'''
    
    with open(task_dir / "grade.py", "w") as f:
        f.write(grade_script)
    
    print(f"✅ Created OpenProblems batch_integration task at {task_dir}")
    print(f"📊 Leaderboard has {len(simple_df)} methods")
    print(f"📈 Full results have {len(leaderboard_df)} method-metric combinations")

if __name__ == "__main__":
    create_batch_integration_task()
```

### Step 4: Manual Execution

1. **Install dependencies**:
```bash
pip install boto3 pyyaml pandas
```

2. **Run the script**:
```bash
python scripts/create_openproblems_poc.py
```

3. **Verify the task**:
```bash
ls biomlbench/tasks/openproblems/batch_integration/
# Should see: config.yaml, leaderboard.csv, description.md, prepare.py, grade.py, results_data.yaml
```

4. **Test with biomlbench**:
```bash
# This should list the new task
python -c "from biomlbench.registry import registry; print(registry.list_task_ids())"
```

## Expected Output

The task will:
- ✅ Display in biomlbench task listing
- ✅ Show a leaderboard with method performance
- ✅ Provide task description and metadata
- ✅ Link back to OpenProblems for full evaluation
- ❌ Cannot be executed by agents (by design)

## Limitations of This PoC

- **Single task only**: Just batch_integration
- **Simplified leaderboard**: Takes first score per method
- **No metric aggregation**: Avoids the scientific validity issues
- **Reference-only**: Cannot be executed, just displays results
- **Hardcoded run**: Uses specific S3 path, not dynamic discovery
- **Minimal error handling**: Focuses on simplicity

## Next Steps If PoC Works

1. Add more OpenProblems tasks
2. Improve metric handling and display
3. Add proper error handling
4. Create automated ingestion pipeline
5. Add validation and testing

This PoC proves the concept works without the complexity and scientific validity issues of the original approach. 