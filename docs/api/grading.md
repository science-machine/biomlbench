# Grading API

The grading module provides evaluation and scoring functionality for agent submissions on biomedical tasks.

## Core Functions

::: biomlbench.grade.grade_jsonl

::: biomlbench.grade.grade_csv

::: biomlbench.grade.validate_submission

::: biomlbench.grade.aggregate_reports

## Grader Classes

::: biomlbench.grade_helpers.Grader

::: biomlbench.grade_helpers.TaskReport

## Single Task Evaluation

```python
from biomlbench.grade import grade_csv
from biomlbench.registry import registry

# Grade a single submission
task = registry.get_task("caco2-wang")
submission_path = Path("submission.csv")

report = grade_csv(submission_path, task)

print(f"Score: {report.score}")
print(f"Medal: {'🥇' if report.gold_medal else '🥈' if report.silver_medal else '🥉' if report.bronze_medal else '❌'}")
print(f"Beats human: {report.beats_human}")
```

## Multi-Task Evaluation

```python
from biomlbench.grade import grade_jsonl
from pathlib import Path

# Grade multiple tasks from submission.jsonl
submission_path = Path("runs/my-run-group/submission.jsonl")
output_dir = Path("results/")

grade_jsonl(submission_path, output_dir)
# Creates timestamped grading report in output_dir
```

## Custom Graders

```python
from biomlbench.grade_helpers import Grader
import pandas as pd

# Define custom metric
def rmse_metric(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """Root Mean Square Error for regression tasks."""
    y_true = answers['target'].values
    y_pred = submission['prediction'].values
    return np.sqrt(np.mean((y_true - y_pred) ** 2))

# Create grader
grader = Grader(name="rmse", grade_fn="my_module.metrics:rmse_metric")

# Use grader
score = grader(submission_df, answers_df)
```

## Biomedical Metrics

BioML-bench includes domain-specific evaluation metrics:

### Medical Imaging
- **AUC-ROC** - Area under ROC curve
- **Precision/Recall** - Binary classification metrics
- **Dice Coefficient** - Segmentation overlap

### Drug Discovery
- **RMSE** - Root mean square error
- **R² Score** - Coefficient of determination  
- **Spearman Correlation** - Rank correlation

### Protein Engineering
- **RMSD** - Root mean square deviation
- **TM-score** - Template modeling score
- **GDT-TS** - Global distance test

## Human Performance Comparison

```python
from biomlbench.grade import calculate_human_performance_metrics

# Compare agent to human baselines
agent_score = 0.89
human_df = pd.read_csv("human_baselines.csv")
is_lower_better = False

beats_human, percentile = calculate_human_performance_metrics(
    agent_score, human_df, is_lower_better
)

print(f"Agent beats human median: {beats_human}")
print(f"Agent percentile vs humans: {percentile}%")
```

## Medal System

BioML-bench uses a Kaggle-style medal system that varies based on leaderboard size:

**For small leaderboards (1-99 teams):**
- **🥇 Gold**: Top 10% of submissions
- **🥈 Silver**: Top 20% (but not gold)  
- **🥉 Bronze**: Top 40% (but not silver/gold)

**For medium leaderboards (100-249 teams):**
- **🥇 Gold**: Top 10 positions (fixed)
- **🥈 Silver**: Top 20% (but not gold)
- **🥉 Bronze**: Top 40% (but not silver/gold)

**For large leaderboards (250-999 teams):**
- **🥇 Gold**: Top (10 + 0.2% of teams) positions
- **🥈 Silver**: Top 50 positions (fixed)
- **🥉 Bronze**: Top 100 positions (fixed)

**For very large leaderboards (1000+ teams):**
- **🥇 Gold**: Top (10 + 0.2% of teams) positions
- **🥈 Silver**: Top 5% of submissions
- **🥉 Bronze**: Top 10% of submissions

Medal thresholds follow the official Kaggle competition progression system.

## Submission Validation

```python
from biomlbench.grade import validate_submission

# Check submission format
is_valid, message = validate_submission(submission_path, task)

if not is_valid:
    print(f"Submission error: {message}")
else:
    print("Submission is valid")
```

## Report Generation

TaskReport objects contain comprehensive evaluation results:

```python
report = grade_csv(submission_path, task)

# Access results
print(f"Task: {report.task_id}")
print(f"Score: {report.score}")
print(f"Gold medal: {report.gold_medal}")
print(f"Above median: {report.above_median}")
print(f"Beats human: {report.beats_human}")
print(f"Human percentile: {report.human_percentile}")

# Convert to dictionary for JSON export
report_dict = report.to_dict()
``` 