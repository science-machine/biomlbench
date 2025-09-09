# Basic Usage

This guide covers the basic workflow for using BioML-bench to evaluate agents on biomedical tasks.

## Overview

The typical BioML-bench workflow involves:

1. **Prepare** task datasets
2. **Run** agents on tasks  
3. **Grade** agent submissions
4. **Analyze** results

## Task Preparation

Before running agents, you need to prepare task datasets:

### Prepare Individual Tasks

```bash
# Prepare a specific task
biomlbench prepare -t polarishub/tdcommons-caco2-wang

# Prepare with options
biomlbench prepare -t manual/histopathologic-cancer-detection --keep-raw
```

### Prepare Multiple Tasks / filtered tasks

**Not implemented yet.**


## Running Agents

### Built-in Agents

BioML-bench includes several reference agents:

```bash
# Dummy agent (for testing)
biomlbench run-agent --agent dummy --task-id polarishub/tdcommons-caco2-wang

# AIDE agent (requires OpenAI API key)
biomlbench run-agent --agent aide --task-id polarishub/tdcommons-caco2-wang
```


## Understanding Agent Outputs

When an agent run completes, BioML-bench creates:

```
runs/
└── 2024-01-15T10-30-00-GMT_run-group_dummy/
    ├── metadata.json              # Run summary
    ├── submission.jsonl           # Auto-generated submission file
    └── polarishub/
        └── <task_id>_abc123/  # Individual task run
            ├── submission/
            │   └── submission.<ext>        # Agent predictions (e.g., `submission.csv`, `submission.h5ad`)
            ├── logs/
            │   └── run.log              # Execution logs
            └── code/                    # Agent code (if available)
```

### Key Files

- **`metadata.json`**: Run metadata and success status
- **`submission.jsonl`**: Automatically generated submission index
- **`submission.<ext>`**: Agent's predictions for each task (e.g., `submission.csv`, `submission.h5ad`)
- **`run.log`**: Detailed execution logs

## Grading and Evaluation

### Grade Agent Runs

The submission file is automatically generated after agent execution:

```bash
# Grade using auto-generated submission file
biomlbench grade \
    --submission runs/2024-01-15T10-30-00-GMT_run-group_dummy/submission.jsonl \
    --output-dir results/
```

### Grade Individual Submissions

For testing or external submissions:

```bash
# Grade a single CSV file
biomlbench grade-sample submission.csv polarishub/tdcommons-caco2-wang

# Grade a single H5AD file
biomlbench grade-sample submission.h5ad openproblems/cell_cell_communication

```



## Working with Results

### Understanding Scores

Each task uses task/domain-specific metrics.

E.g., Caco2-Wang uses MAE, while Histopathologic Cancer Detection uses AUC-ROC.

### Medal System

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

This follows the official Kaggle competition progression system.

### Human Baseline Comparison

Many tasks include human expert performance for context:

```json
{
  "task_id": "manual/histopathologic-cancer-detection",
  "score": 0.89,
  "beats_human": true,
  "human_percentile": 85.5
}
```

## Advanced Usage

### Custom Container Configuration

Override Docker settings for specific requirements:

```json
// custom_config.json
{
    "mem_limit": "8g",
    "shm_size": "4g", 
    "nano_cpus": 4000000000,
    "gpus": -1
}
```

```bash
biomlbench run-agent \
    --agent my-agent \
    --task-id polarishub/tdcommons-caco2-wang \
    --container-config custom_config.json
```


### Task Development Workflow

For full details, see [developer/adding_tasks.md](../docs/developer/adding_tasks.md)

```bash
# 1. Create new task structure
mkdir -p biomlbench/tasks/my-source/my-new-task

# 2. Test task preparation
biomlbench prepare -t my-source/my-new-task

# 3. Test with dummy agent
biomlbench run-agent --agent dummy --task-id my-source/my-new-task

# 4. Grade test submission (if CSV format)
biomlbench grade-sample runs/dummy/my-source/my-new-task/submission.csv my-source/my-new-task
# Grade test submission (if H5AD format)
biomlbench grade-sample runs/dummy/my-source/my-new-task/submission.h5ad my-source/my-new-task
```

## CLI Reference

### Main Commands

- **`prepare`**: Download and prepare task datasets
- **`run-agent`**: Execute agents on tasks
- **`grade`**: Evaluate submissions  
- **`grade-sample`**: Test individual submissions
- **`run-baseline`**: Generate baseline comparisons

### Getting Help

```bash
# Main help
biomlbench --help

# Command-specific help
biomlbench prepare --help
biomlbench run-agent --help
biomlbench grade --help
```

## Common Workflows

### Quick Testing

```bash
# Fast workflow for testing
biomlbench prepare -t polarishub/tdcommons-caco2-wang
biomlbench run-agent --agent dummy --task-id polarishub/tdcommons-caco2-wang
biomlbench grade --submission runs/*/submission.jsonl --output-dir results/
```

### Full Evaluation

```bash
# Comprehensive evaluation (when more tasks are available)
# biomlbench prepare --all
# biomlbench run-agent --agent my-agent --task-list experiments/splits/all.txt --n-workers 4
# biomlbench grade --submission runs/*/submission.jsonl --output-dir results/

# Current example with available tasks
biomlbench prepare -t polarishub/tdcommons-caco2-wang -t manual/histopathologic-cancer-detection
biomlbench run-agent --agent my-agent --task-id polarishub/tdcommons-caco2-wang
biomlbench grade --submission runs/*/submission.jsonl --output-dir results/
```
