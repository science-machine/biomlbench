# Grading API

Evaluation and scoring functionality for agent submissions on biomedical tasks.

::: biomlbench.grade
    options:
      show_source: false
      heading_level: 2

::: biomlbench.grade_helpers
    options:
      show_source: false
      heading_level: 2

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

## Usage Examples

### Single Task Evaluation

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

### Multi-Task Evaluation

```python
from biomlbench.grade import grade_jsonl

# Grade multiple tasks from submission.jsonl
grade_jsonl(
    path_to_submissions=Path("runs/my-run-group/submission.jsonl"),
    output_dir=Path("results/")
)
``` 