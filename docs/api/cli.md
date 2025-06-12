# CLI Reference

BioML-bench's command-line interface provides comprehensive tools for managing biomedical benchmark tasks and running agent evaluations.

## Main Entry Point

::: biomlbench.cli.main

## Commands Overview

The CLI provides several subcommands for different operations:

### Task Management
- **`prepare`** - Download and prepare task datasets
- **`dev download-leaderboard`** - Download competition leaderboards
- **`dev prepare-human-baselines`** - Extract human performance baselines

### Agent Execution  
- **`run-agent`** - Execute AI agents on biomedical tasks
- **`run-baseline`** - Run baseline algorithms for comparison

### Evaluation
- **`grade`** - Evaluate agent submissions across multiple tasks
- **`grade-sample`** - Test individual task submissions

## Command Examples

### Prepare Tasks

```bash
# Single task
biomlbench prepare -t caco2-wang

# Multiple tasks by domain
biomlbench prepare --domain oncology

# All easy/medium tasks
biomlbench prepare --lite

# Custom data directory
biomlbench prepare -t caco2-wang --data-dir /custom/path
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

# With custom container configuration
biomlbench run-agent \
    --agent my-agent \
    --task-id caco2-wang \
    --container-config custom_docker.json \
    --retain-container
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

### Run Baselines

```bash
# Specific baseline
biomlbench run-baseline caco2-wang --baseline linear

# All baselines
biomlbench run-baseline caco2-wang --baseline all

# Custom output directory
biomlbench run-baseline caco2-wang --baseline rf --output-dir baselines/
```

## Command Line Arguments

### Global Options

All commands support these options:

- **`--help`** - Show command help
- **`--data-dir`** - Custom data directory (overrides default cache location)

### Task Filtering

For `prepare` command:

- **`-t, --task-id`** - Specific task ID
- **`--all`** - All available tasks  
- **`--lite`** - Easy/medium difficulty tasks only
- **`--domain`** - Tasks from specific biomedical domain
- **`--task-type`** - Tasks of specific type (medical_imaging, drug_discovery, etc.)
- **`-l, --list`** - Tasks from file (one per line)

### Agent Execution

For `run-agent` command:

- **`--agent`** - Agent ID (required)
- **`--task-id`** - Single task to run
- **`--task-list`** - File with multiple tasks
- **`--n-workers`** - Parallel workers (default: 1)
- **`--n-seeds`** - Random seeds per task (default: 1)
- **`--container-config`** - Custom Docker configuration JSON
- **`--retain-container`** - Keep containers after run
- **`--output-dir`** - Custom output directory

### Evaluation Options

For `grade` command:

- **`--submission`** - Path to submission JSONL file (required)
- **`--output-dir`** - Results directory (required)

For `grade-sample` command:

- **`submission`** - Path to CSV file (positional)
- **`task_id`** - Task ID (positional)

## Environment Variables

BioML-bench respects these environment variables:

### Agent Configuration  
- **`OPENAI_API_KEY`** - For AIDE agent
- **`I_ACCEPT_RUNNING_PRIVILEGED_CONTAINERS`** - Allow privileged containers (set to "true")

### Kaggle API
Kaggle authentication uses the standard Kaggle API configuration file (`~/.kaggle/kaggle.json`). See the [Kaggle API documentation](https://github.com/Kaggle/kaggle-api) for setup instructions.

## Configuration Files

### Container Configuration

Custom Docker settings (JSON format):

```json
{
    "mem_limit": "8g",
    "shm_size": "4g",
    "nano_cpus": 4000000000,
    "gpus": -1,
    "runtime": "sysbox-runc"
}
```

### Task Lists

Task list files contain one task ID per line:

```text
caco2-wang
histopathologic-cancer-detection
# Comments are supported
```

## Exit Codes

BioML-bench uses standard exit codes:

- **`0`** - Success
- **`1`** - General error
- **`2`** - Invalid command line arguments
- **`130`** - Interrupted by user (Ctrl+C)

## Logging

CLI output includes structured logging:

```
[2024-01-15 10:30:00,123] [cli.py:45] Running agent 'dummy' on tasks: ['caco2-wang']
[2024-01-15 10:30:05,456] [agents.py:78] Launching run group: 2024-01-15T10-30-00-GMT_run-group_dummy
[2024-01-15 10:32:15,789] [utils.py:234] Generated submission file: runs/.../submission.jsonl
```

Log levels can be controlled via Python logging configuration.

## Shell Completion

Generate shell completion for improved CLI experience:

```bash
# Bash
biomlbench --help | grep -A 20 "completion"

# Zsh  
compinit && complete -o default -F _biomlbench biomlbench

# Fish
biomlbench --help | fish_completion
``` 