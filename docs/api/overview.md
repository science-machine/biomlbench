# API Reference Overview

This section provides comprehensive API documentation for BioML-bench's Python modules and classes.

## Core Modules

### [CLI](cli.md) - `biomlbench.cli`
Command-line interface and argument parsing for all BioML-bench commands.

### [Registry](registry.md) - `biomlbench.registry`  
Task registry system for discovering, loading, and managing biomedical benchmark tasks.

### [Tasks](tasks.md) - `biomlbench.data`
Task preparation, data downloading, and dataset management functionality.

### [Agents](agents.md) - `biomlbench.agents`
Agent execution system for running AI agents in containerized environments.

### [Grading](grading.md) - `biomlbench.grade`
Evaluation and scoring system for agent submissions with biomedical metrics.

### [Data Sources](data_sources.md) - `biomlbench.data_sources`
Pluggable data source system supporting Kaggle, Polaris, and custom sources.

### [Utilities](utils.md) - `biomlbench.utils`
Shared utility functions for file operations, logging, and data processing.

## Key Classes

### Task
Represents a biomedical benchmark task with metadata, data paths, and evaluation configuration.

```python
from biomlbench import registry

task = registry.get_task("caco2-wang")
print(f"Task: {task.name}")
print(f"Domain: {task.domain}")
print(f"Type: {task.task_type}")
```

### Registry
Central registry for discovering and managing tasks across different biomedical domains.

```python
from biomlbench.registry import registry

# List all available tasks
tasks = registry.list_task_ids()

# Get tasks by domain (not implemented yet)
# oncology_tasks = registry.get_tasks_by_domain("oncology")

# Get tasks by type (not implemented yet)
# imaging_tasks = registry.get_tasks_by_type("medical_imaging")
```

### Agent  
Represents an AI agent with execution configuration and container settings.

```python
from agents.registry import registry as agent_registry

agent = agent_registry.get_agent("dummy")
print(f"Agent: {agent.name}")
print(f"Privileged: {agent.privileged}")
```

### Grader
Handles task-specific evaluation metrics and scoring logic. 

**Note:** This works but we don't really use it anywhere yet. Possibly remove unless we really have some more complex metrics needed across tasks. (TODO)

```python
from biomlbench.grade_helpers import Grader

grader = Grader(name="auc-roc", grade_fn="biomlbench.metrics:average_precision_at_k")
score = grader(submission_df, answers_df)
```

## Data Source Architecture

BioML-bench uses a pluggable data source system:

```python
from biomlbench.data_sources import create_data_source, list_available_sources
from pathlib import Path

# List available data source types
available_sources = list_available_sources()
print(f"Available sources: {available_sources}")

# Create Kaggle data source
kaggle_source = create_data_source("kaggle")

# Download task data
kaggle_source.download(
    source_config={"competition_id": "histopathologic-cancer-detection"},
    data_dir=Path("test/")
)

# Get leaderboard
leaderboard = kaggle_source.get_leaderboard(
    source_config={"competition_id": "histopathologic-cancer-detection"}
)
```

## Task Configuration

Tasks are configured via YAML files with biomedical-specific metadata:

```yaml
id: caco2-wang
name: "Caco-2 Permeability Prediction"
task_type: drug_discovery
domain: pharmacokinetics
difficulty: medium

dataset:
  answers: caco2-wang/prepared/private/answers.csv
  sample_submission: caco2-wang/prepared/public/sample_submission.csv

grader:
  name: rmse
  grade_fn: biomlbench.tasks.caco2-wang.grade:grade

biomedical_metadata:
  modality: "molecular_properties"
  data_type: "regression"
  clinical_relevance: "drug_absorption"
```

## Error Handling

BioML-bench defines custom exceptions for different error conditions:

```python
from biomlbench.data_sources.base import DataSourceError
from biomlbench.grade_helpers import InvalidSubmissionError

try:
    score = grader(submission, answers)
except InvalidSubmissionError as e:
    logger.warning(f"Invalid submission: {e}")
except DataSourceError as e:
    logger.error(f"Data source error: {e}")
```

## Logging

BioML-bench uses structured logging throughout:

```python
from biomlbench.utils import get_logger

logger = get_logger(__name__)
logger.info("Task preparation starting")
logger.warning("Missing optional dependency")
logger.error("Task preparation failed")
```

## Extension Points

### Adding Data Sources

Implement the `DataSource` interface:

```python
from biomlbench.data_sources.base import DataSource

class MyDataSource(DataSource):
    def download(self, source_config, data_dir):
        # Implementation
        pass
    
    def get_leaderboard(self, source_config):
        # Implementation  
        pass
```

See examples for kaggle: `biomlbench/data_sources/kaggle.py` and polaris: `biomlbench/data_sources/polaris.py`.

### Adding Metrics

Implement grading functions:

```python
def my_metric(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """Custom biomedical metric."""
    # Implementation
    return score
```

### Adding Tasks

See the documentation on [adding tasks](../../docs/developer/adding_tasks.md).

### Adding Agents

See the documentation on [creating agents](../../docs/developer/creating_agents.md).

## Type Hints

BioML-bench uses type hints throughout for better IDE support:

```python
from typing import List, Optional, Dict, Any
from pathlib import Path
import pandas as pd

def process_task(
    task_id: str,
    data_dir: Optional[Path] = None
) -> Dict[str, Any]:
    """Process a biomedical task with type safety."""
    pass
``` 