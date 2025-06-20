# Adding Tasks

Guide for adding new biomedical benchmark tasks to BioML-bench.

## Task Requirements

Every task needs:

- Pre-split data for training and testing
- Clear evaluation metrics
- A description of the task

**Note:** If adding a benchmarks from a new database, you'll also need to add a new data source module in `biomlbench/data_sources/`. See examples in `biomlbench/data_sources/kaggle.py` and `biomlbench/data_sources/polaris.py`.

## Implementation Steps

1. **Create task directory structure**
2. **Configure task metadata**
3. **Implement data preparation**
4. **Define evaluation logic**
5. **Write task description**
6. **Test and validate**

## Structure of task directory

```
tasks/
├── my-source/                   # Data source folder (e.g., polarishub, manual)
│   └── my-biomedical-task/      # Task directory
│       ├── config.yaml          # Task configuration
│       ├── description.md       # Task description
│       ├── prepare.py           # Data preparation script
│       ├── grade.py             # Evaluation logic
│       └── prepared/            # Generated during preparation
│           ├── public/          # Public data
│           │   ├── train.csv   # Training data
│           │   ├── test_features.csv # Test features
│           │   └── sample_submission.csv # Example submission
│           └── private/         # Private data
│               └── answers.csv  # Test set answers

```

## Task Configuration (`config.yaml`)

```yaml
id: my-source/my-biomedical-task
name: "My Biomedical Task"
task_type: drug_discovery  # or medical_imaging, protein_engineering
domain: pharmacokinetics   # specific biomedical domain
difficulty: medium         # easy, medium, hard

data_source:
  type: kaggle            # or polaris, custom
  competition_id: my-task

dataset:
  answers: my-source/my-biomedical-task/prepared/private/answers.csv
  sample_submission: my-source/my-biomedical-task/prepared/public/sample_submission.csv

grader:
  name: rmse
  grade_fn: biomlbench.tasks.my-source.my-biomedical-task.grade:grade

preparer: biomlbench.tasks.my-source.my-biomedical-task.prepare:prepare

biomedical_metadata:
  modality: "molecular_properties"
  organ_system: "liver"
  data_type: "regression"
  clinical_relevance: "drug_metabolism"
```

## Data Preparation (`prepare.py`)

See example `biomlbench/tasks/polarishub/tdcommons-caco2-wang/prepare.py`.

```python
from pathlib import Path
import pandas as pd

def prepare(task_dir: Path, raw_dir: Path, public_dir: Path, private_dir: Path) -> Path:
    """Prepare task data with train/test splits."""
    
    # Download and process raw data
    # Create train.csv, test_features.csv, sample_submission.csv
    # Generate private answers.csv
    
    return public_dir
```

## Evaluation Logic (`grade.py`)

See example `biomlbench/tasks/polarishub/tdcommons-caco2-wang/grade.py`.

```python
import pandas as pd
import numpy as np

def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """Calculate task-specific metric."""
    
    y_true = answers['label'].values
    y_pred = submission['label'].values
    
    # Implement domain-specific metric
    return np.sqrt(np.mean((y_true - y_pred) ** 2))
```

## Testing New Tasks

```bash
# Test preparation
biomlbench prepare -t my-source/my-task

# Test with dummy agent
biomlbench run-agent --agent dummy --task-id my-source/my-task

# Validate submission
biomlbench grade --submission /path/to/submission.jsonl --output-dir /path/to/output/dir
``` 