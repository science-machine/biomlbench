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
biomlbench prepare -t caco2-wang

# Prepare with options
biomlbench prepare -t histopathologic-cancer-detection --keep-raw
```

### Prepare Multiple Tasks

```bash
# Prepare by difficulty
biomlbench prepare --lite  # Easy/medium tasks only

# Prepare by domain
biomlbench prepare --domain oncology

# Prepare by task type  
biomlbench prepare --task-type medical_imaging

# Prepare from a list
biomlbench prepare --list experiments/splits/drug-discovery.txt

# Prepare all tasks
biomlbench prepare --all
```

### Available Task Filters

- **`--lite`**: Low difficulty tasks for quick testing
- **`--domain`**: Filter by biomedical domain (oncology, cardiology, etc.)
- **`--task-type`**: Filter by task type (medical_imaging, drug_discovery, etc.)
- **`--list`**: Use a custom task list file

## Running Agents

### Built-in Agents

BioML-bench includes several reference agents:

```bash
# Dummy agent (for testing)
biomlbench run-agent --agent dummy --task-id caco2-wang

# AIDE agent (requires OpenAI API key)
biomlbench run-agent --agent aide --task-id caco2-wang
```

### Single Task Execution

```bash
# Run agent on one task
biomlbench run-agent --agent dummy --task-id caco2-wang

# With custom settings
biomlbench run-agent \
    --agent dummy \
    --task-id caco2-wang \
    --n-seeds 3 \
    --retain-container
```

### Multi-Task Execution

```bash
# Run on multiple tasks from a list
biomlbench run-agent --agent dummy --task-list experiments/splits/all.txt

# Parallel execution
biomlbench run-agent \
    --agent dummy \
    --task-list experiments/splits/medical-imaging.txt \
    --n-workers 4 \
    --n-seeds 2
```

### Agent Execution Options

- **`--n-workers`**: Number of parallel workers (default: 1)
- **`--n-seeds`**: Random seeds per task (default: 1)  
- **`--retain-container`**: Keep containers for debugging
- **`--container-config`**: Custom Docker configuration
- **`--data-dir`**: Custom data directory

## Understanding Agent Outputs

When an agent run completes, BioML-bench creates:

```
runs/
└── 2024-01-15T10-30-00-GMT_run-group_dummy/
    ├── metadata.json              # Run summary
    ├── submission.jsonl           # Auto-generated submission file
    ├── caco2-wang_abc123/         # Individual task run
    │   ├── submission/
    │   │   └── submission.csv     # Agent predictions
    │   ├── logs/
    │   │   └── run.log           # Execution logs
    │   └── code/                 # Agent code (if available)
    └── task2_def456/             # Additional task runs...
```

### Key Files

- **`metadata.json`**: Run metadata and success status
- **`submission.jsonl`**: Automatically generated submission index
- **`submission.csv`**: Agent's predictions for each task
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
biomlbench grade-sample submission.csv caco2-wang
```

### Baseline Comparisons

Run baselines to establish performance benchmarks:

```bash
# Run specific baseline
biomlbench run-baseline caco2-wang --baseline linear

# Run all available baselines
biomlbench run-baseline caco2-wang --baseline all

# Grade baseline results
biomlbench grade --submission baseline_submissions/submission.jsonl --output-dir results/
```

## Working with Results

### Understanding Scores

Each task uses domain-specific metrics:

- **Medical Imaging**: AUC-ROC, Precision, Recall
- **Drug Discovery**: RMSE, R², Spearman correlation  
- **Protein Engineering**: RMSD, TM-score

### Medal System

BioML-bench uses a Kaggle-style medal system based on leaderboard percentiles:

- **🥇 Gold**: Top 10% or top 10 submissions
- **🥈 Silver**: Top 20% (but not gold)
- **🥉 Bronze**: Top 40% (but not silver/gold)

### Human Baseline Comparison

Many tasks include human expert performance for context:

```json
{
  "task_id": "histopathologic-cancer-detection",
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
    --task-id caco2-wang \
    --container-config custom_config.json
```

### Custom Data Directory

Use a different data storage location:

```bash
biomlbench prepare -t caco2-wang --data-dir /custom/data/path

biomlbench run-agent \
    --agent dummy \
    --task-id caco2-wang \
    --data-dir /custom/data/path
```

### Task Development Workflow

```bash
# 1. Create new task structure
mkdir -p biomlbench/tasks/my-new-task

# 2. Test task preparation
biomlbench prepare -t my-new-task

# 3. Test with dummy agent
biomlbench run-agent --agent dummy --task-id my-new-task

# 4. Grade test submission
biomlbench grade-sample runs/dummy/my-new-task/submission.csv my-new-task
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
biomlbench prepare -t caco2-wang
biomlbench run-agent --agent dummy --task-id caco2-wang
biomlbench grade --submission runs/*/submission.jsonl --output-dir results/
```

### Full Evaluation

```bash
# Comprehensive evaluation
biomlbench prepare --all
biomlbench run-agent --agent my-agent --task-list experiments/splits/all.txt --n-workers 4
biomlbench grade --submission runs/*/submission.jsonl --output-dir results/
```

### Development Cycle

```bash
# Iterative development
biomlbench prepare -t my-task
biomlbench run-agent --agent my-agent --task-id my-task --retain-container
# Debug containers, fix issues, repeat
``` 