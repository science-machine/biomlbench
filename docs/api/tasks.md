# Tasks and Data API

The data module handles task preparation, dataset downloading, and data management for biomedical benchmark tasks.

## Core Functions

::: biomlbench.data.download_and_prepare_dataset

::: biomlbench.data.is_dataset_prepared

::: biomlbench.data.ensure_leaderboard_exists

::: biomlbench.data.get_leaderboard

::: biomlbench.data.prepare_human_baselines

## Task Preparation Workflow

```python
from biomlbench.registry import registry
from biomlbench.data import download_and_prepare_dataset, is_dataset_prepared

# Get task from registry
task = registry.get_task("caco2-wang")

# Check if already prepared
if not is_dataset_prepared(task):
    # Download and prepare dataset
    download_and_prepare_dataset(
        task=task,
        keep_raw=False,
        skip_verification=False
    )

print(f"Task prepared: {task.name}")
```

## Data Directory Structure

After preparation, tasks have the following structure:

```
data/
└── task-id/
    ├── raw/                    # Raw downloaded data
    ├── prepared/
    │   ├── public/            # Data accessible to agents
    │   │   ├── train.csv
    │   │   ├── test_features.csv
    │   │   ├── sample_submission.csv
    │   │   └── description.md
    │   └── private/           # Private evaluation data
    │       └── answers.csv
    └── checksums.yaml         # Data integrity verification
```

## Data Source Integration

Tasks can use different data sources:

```python
# Kaggle competition data
task_config = {
    "data_source": {
        "type": "kaggle",
        "competition_id": "histopathologic-cancer-detection"
    }
}

# Polaris molecular datasets
task_config = {
    "data_source": {
        "type": "polaris",
        "dataset_id": "Caco2_Wang"
    }
}
```

## Human Baseline Management

```python
from biomlbench.data import prepare_human_baselines

# Extract human performance data
task = registry.get_task("histopathologic-cancer-detection")
prepare_human_baselines(task, force=True)

# Human baselines are saved as CSV files
human_baselines_path = task.public_dir / "human_baselines.csv"
```

## Verification and Validation

The module provides data integrity checking:

- **Checksums** - Verify downloaded data integrity
- **Format validation** - Ensure CSV files have correct structure
- **Completeness checks** - Verify all required files are present 