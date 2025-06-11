# Adding Tasks

Detailed guide for adding new biomedical benchmark tasks to BioML-bench.

## Task Requirements

Every task needs:
- Biomedical relevance and scientific validity
- Clear evaluation metrics
- Sufficient data for training and testing
- Proper train/test splits to prevent data leakage

## Implementation Steps

1. **Create task directory structure**
2. **Configure task metadata**
3. **Implement data preparation**
4. **Define evaluation logic**
5. **Write task description**
6. **Test and validate**

## Task Configuration (`config.yaml`)

```yaml
id: my-biomedical-task
name: "My Biomedical Task"
task_type: drug_discovery  # or medical_imaging, protein_engineering
domain: pharmacokinetics   # specific biomedical domain
difficulty: medium         # easy, medium, hard

data_source:
  type: kaggle            # or polaris, custom
  competition_id: my-task

dataset:
  answers: my-task/prepared/private/answers.csv
  sample_submission: my-task/prepared/public/sample_submission.csv

grader:
  name: rmse
  grade_fn: biomlbench.tasks.my-task.grade:grade

preparer: biomlbench.tasks.my-task.prepare:prepare

biomedical_metadata:
  modality: "molecular_properties"
  organ_system: "liver"
  data_type: "regression"
  clinical_relevance: "drug_metabolism"
```

## Data Preparation (`prepare.py`)

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

```python
import pandas as pd
import numpy as np

def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """Calculate task-specific metric."""
    
    y_true = answers['target'].values
    y_pred = submission['prediction'].values
    
    # Implement domain-specific metric
    return np.sqrt(np.mean((y_true - y_pred) ** 2))
```

## Testing New Tasks

```bash
# Test preparation
biomlbench prepare -t my-task

# Test with dummy agent
biomlbench run-agent --agent dummy --task-id my-task

# Validate submission
biomlbench grade-sample submission.csv my-task
``` 