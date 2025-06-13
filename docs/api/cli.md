# CLI Reference

BioML-bench's command-line interface for managing biomedical benchmark tasks and running agent evaluations.

::: biomlbench.cli
    options:
      show_source: false
      heading_level: 2

## Environment Variables

BioML-bench respects these environment variables:

### Agent Configuration  
- **`OPENAI_API_KEY`** - For AIDE agent
- **`I_ACCEPT_RUNNING_PRIVILEGED_CONTAINERS`** - Allow privileged containers (set to "true")

### Kaggle API
Kaggle authentication uses the standard Kaggle API configuration file (`~/.kaggle/kaggle.json`). See the [Kaggle API documentation](https://github.com/Kaggle/kaggle-api) for setup instructions.

## Usage Examples

### Prepare Tasks

```bash
# Single task
biomlbench prepare -t caco2-wang

# Multiple tasks by domain
biomlbench prepare --domain admet

# All easy/medium tasks
biomlbench prepare --lite
```

### Run Agents

```bash
# Single task
biomlbench run-agent --agent dummy --task-id caco2-wang

# Multiple tasks with parallelization (example with available split)
biomlbench run-agent \
    --agent aide \
    --task-list experiments/splits/caco2-wang.txt \
    --n-workers 4 \
    --n-seeds 2
```

### Grade Submissions

```bash
# Grade agent run (submission.jsonl auto-generated)
biomlbench grade \
    --submission runs/2024-01-15_run-group_dummy/submission.jsonl \
    --output-dir results/

# Test single submission
biomlbench grade-sample submission.csv caco2-wang
``` 