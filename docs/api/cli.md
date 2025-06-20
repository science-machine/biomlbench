# CLI Reference

BioML-bench's command-line interface for managing biomedical benchmark tasks and running agent evaluations.

::: biomlbench.cli
    options:
      show_source: false
      heading_level: 2

## Main Help

```
usage: cli.py [-h] {prepare,grade,grade-sample,dev,run-baseline,run-agent} ...

Runs agents on biomedical ML tasks.

positional arguments:
  {prepare,grade,grade-sample,dev,run-baseline,run-agent}
                        Sub-command to run.
    prepare             Download and prepare tasks for the BioML-bench
                        dataset.
    grade               Grade a submission to the eval, comprising of several
                        task submissions
    grade-sample        Grade a single sample (task) in the eval
    dev                 Developer tools for extending BioML-bench.
    run-baseline        Run a baseline agent on a biomedical task
    run-agent           Run an AI agent on biomedical tasks

options:
  -h, --help            show this help message and exit

```

## Commands

### `biomlbench prepare`

```
usage: cli.py prepare [-h] [-t TASK_ID] [-a] [--lite] [-l LIST]
                      [--domain DOMAIN] [--task-type TASK_TYPE] [--keep-raw]
                      [--data-dir DATA_DIR] [--overwrite-checksums]
                      [--overwrite-leaderboard] [--skip-verification]

options:
  -h, --help            show this help message and exit
  -t TASK_ID, --task-id TASK_ID
                        ID of the task to prepare. Valid options:
                        ['polarishub/tdcommons-caco2-wang', 'manual/histopathologic-cancer-detection']
  -a, --all             Prepare all tasks.
  --lite                Prepare all the low difficulty tasks (BioML-bench
                        Lite).
  -l LIST, --list LIST  Prepare a list of tasks specified line by line in a
                        text file.
  --domain DOMAIN       Prepare all tasks for a specific biomedical domain
                        (e.g., oncology, drug_discovery).
  --task-type TASK_TYPE
                        Prepare all tasks of a specific type (e.g.,
                        medical_imaging, protein_engineering).
  --keep-raw            Keep the raw task files after the task has been
                        prepared.
  --data-dir DATA_DIR   Path to the directory where the data will be stored.
  --overwrite-checksums
                        [For Developers] Overwrite the checksums file for the
                        task.
  --overwrite-leaderboard
                        [For Developers] Overwrite the leaderboard file for
                        the task.
  --skip-verification   [For Developers] Skip the verification of the
                        checksums.

```

### `biomlbench run-agent`

```
usage: cli.py run-agent [-h] --agent AGENT [--task-id TASK_ID]
                        [--task-list TASK_LIST] [--n-workers N_WORKERS]
                        [--n-seeds N_SEEDS]
                        [--container-config CONTAINER_CONFIG]
                        [--retain-container] [--data-dir DATA_DIR]
                        [--output-dir OUTPUT_DIR]

options:
  -h, --help            show this help message and exit
  --agent AGENT         Agent ID to run (e.g., dummy, aide, aide/dev)
  --task-id TASK_ID     Single task ID to run. Valid options: ['polarishub/tdcommons-caco2-wang',
                        'manual/histopathologic-cancer-detection']
  --task-list TASK_LIST
                        Path to text file with task IDs (one per line) for
                        multi-task runs
  --n-workers N_WORKERS
                        Number of parallel workers for multi-task runs
  --n-seeds N_SEEDS     Number of random seeds to run per task
  --container-config CONTAINER_CONFIG
                        Path to JSON file with Docker container configuration
  --retain-container    Keep container after run for debugging
  --data-dir DATA_DIR   Path to the directory where task data is stored
  --output-dir OUTPUT_DIR
                        Directory to save agent run outputs

```

### `biomlbench grade`

```
usage: cli.py grade [-h] --submission SUBMISSION --output-dir OUTPUT_DIR
                    [--data-dir DATA_DIR]

options:
  -h, --help            show this help message and exit
  --submission SUBMISSION
                        Path to the JSONL file of submissions. Refer to
                        README.md#submission-format for the required format.
  --output-dir OUTPUT_DIR
                        Path to the directory where the evaluation metrics
                        will be saved.
  --data-dir DATA_DIR   Path to the directory where the data used for grading
                        is stored.

```

### `biomlbench grade-sample`

```
usage: cli.py grade-sample [-h] [--data-dir DATA_DIR] submission task_id

positional arguments:
  submission           Path to the submission CSV file.
  task_id              ID of the task to grade. Valid options: ['polarishub/tdcommons-caco2-wang',
                       'manual/histopathologic-cancer-detection']

options:
  -h, --help           show this help message and exit
  --data-dir DATA_DIR  Path to the directory where the data will be stored.

```

### `biomlbench run-baseline`

```
usage: cli.py run-baseline [-h] [--baseline BASELINE]
                           [--output-dir OUTPUT_DIR] [--data-dir DATA_DIR]
                           [--seed SEED]
                           task_id

positional arguments:
  task_id               ID of the task to run baseline on. Valid options:
                        ['polarishub/tdcommons-caco2-wang', 'manual/histopathologic-cancer-detection']

options:
  -h, --help            show this help message and exit
  --baseline BASELINE   Baseline to run (simple, random, or task-specific
                        types like linear, rf, fingerprint). Use 'all' to run
                        all available baselines for the task.
  --output-dir OUTPUT_DIR
                        Directory to save baseline submissions
  --data-dir DATA_DIR   Path to the directory where the data is stored.
  --seed SEED           Random seed for reproducible baselines

```

### `biomlbench dev`

```
usage: cli.py dev [-h] {download-leaderboard,prepare-human-baselines} ...

positional arguments:
  {download-leaderboard,prepare-human-baselines}
                        Developer command to run.
    download-leaderboard
                        Download the leaderboard for a task.
    prepare-human-baselines
                        Prepare human baseline data for tasks.

options:
  -h, --help            show this help message and exit

```

#### `biomlbench dev download-leaderboard`

```
usage: cli.py dev download-leaderboard [-h] [-t TASK_ID] [--all] [--force]

options:
  -h, --help            show this help message and exit
  -t TASK_ID, --task-id TASK_ID
                        Name of the task to download the leaderboard for.
                        Valid options: ['polarishub/tdcommons-caco2-wang', 'manual/histopathologic-cancer-
                        detection']
  --all                 Download the leaderboard for all tasks.
  --force               Force download the leaderboard, even if it already
                        exists.

```

#### `biomlbench dev prepare-human-baselines`

```
usage: cli.py dev prepare-human-baselines [-h] [-t TASK_ID] [--all] [--force]

options:
  -h, --help            show this help message and exit
  -t TASK_ID, --task-id TASK_ID
                        Name of the task to prepare human baselines for. Valid
                        options: ['polarishub/tdcommons-caco2-wang', 'manual/histopathologic-cancer-
                        detection']
  --all                 Prepare human baselines for all tasks.
  --force               Force re-extraction of human baselines, even if they
                        already exist.

```

## Usage Examples

### Task Preparation

```bash
# Prepare a specific task
biomlbench prepare -t polarishub/tdcommons-caco2-wang

# Prepare all tasks in a domain
biomlbench prepare --domain admet

# Prepare multiple tasks from a file
biomlbench prepare --list experiments/splits/caco2-wang.txt

# Prepare all low-difficulty tasks
biomlbench prepare --lite
```

### Agent Execution

```bash
# Run agent on single task
biomlbench run-agent --agent dummy --task-id polarishub/tdcommons-caco2-wang

# Run agent on multiple tasks with parallel workers
biomlbench run-agent \
    --agent aide \
    --task-list experiments/splits/caco2-wang.txt \
    --n-workers 4 \
    --n-seeds 3

# Run with custom container configuration
biomlbench run-agent \
    --agent aide \
    --task-id polarishub/tdcommons-caco2-wang \
    --container-config custom_config.json \
    --retain-container
```

### Evaluation and Grading

```bash
# Grade multiple task submissions
biomlbench grade \
    --submission runs/my-run-group/submission.jsonl \
    --output-dir results/

# Grade single task submission
biomlbench grade-sample submission.csv polarishub/tdcommons-caco2-wang

# Run and grade baselines
biomlbench run-baseline polarishub/tdcommons-caco2-wang --baseline all
biomlbench grade \
    --submission baseline_submissions/submission.jsonl \
    --output-dir baseline_results/
```

### Development Commands

```bash
# Download leaderboards for all tasks
biomlbench dev download-leaderboard --all

# Prepare human baselines for a specific task
biomlbench dev prepare-human-baselines -t polarishub/tdcommons-caco2-wang --force
```

## Environment Variables

BioML-bench respects these environment variables:

### Agent Configuration  
- **`OPENAI_API_KEY`** - API key for AIDE agent
- **`I_ACCEPT_RUNNING_PRIVILEGED_CONTAINERS`** - Set to "true" to allow privileged containers

### Kaggle Integration
Kaggle authentication uses the standard configuration file at `~/.kaggle/kaggle.json`. See the [Kaggle API documentation](https://github.com/Kaggle/kaggle-api) for setup instructions.

## Common Patterns

### Full Workflow Example

```bash
# 1. Prepare tasks
biomlbench prepare --lite

# 2. Run agent
biomlbench run-agent --agent dummy --task-list experiments/splits/polarishub-tdcommons-caco2-wang.txt

# 3. Grade results (submission.jsonl is auto-generated)
biomlbench grade \
    --submission runs/latest-run-group/submission.jsonl \
    --output-dir results/
```

### Debugging Agent Issues

```bash
# Run with container retention for debugging
biomlbench run-agent \
    --agent my-agent \
    --task-id polarishub/tdcommons-caco2-wang \
    --retain-container

# Check logs in the run directory
ls runs/latest-run-group/*/logs/
```
