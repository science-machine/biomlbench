# Registry API

The registry system manages discovery, loading, and organization of biomedical benchmark tasks.

## Registry Class

::: biomlbench.registry.Registry

## Task Class

::: biomlbench.registry.Task

## Global Registry Instance

The default registry instance provides immediate access to tasks:

```python
from biomlbench.registry import registry

# Get task information
task = registry.get_task("caco2-wang")
print(f"Task: {task.name} ({task.domain})")

# List all tasks
all_tasks = registry.list_task_ids()
print(f"Available tasks: {len(all_tasks)}")

# Filter by domain
oncology_tasks = registry.get_tasks_by_domain("oncology")
drug_discovery_tasks = registry.get_tasks_by_domain("drug_discovery")

# Filter by task type
imaging_tasks = registry.get_tasks_by_type("medical_imaging")
protein_tasks = registry.get_tasks_by_type("protein_engineering")
```

## Task Discovery

Tasks are automatically discovered from the `biomlbench/tasks/` directory:

```
biomlbench/tasks/
├── caco2-wang/
│   ├── config.yaml          # Task configuration
│   ├── prepare.py           # Data preparation
│   ├── grade.py            # Evaluation logic
│   └── description.md      # Task description
└── histopathologic-cancer-detection/
    ├── config.yaml
    ├── prepare.py
    ├── grade.py
    └── description.md
```

### Task Configuration Format

Each task requires a `config.yaml` file:

```yaml
id: caco2-wang
name: "Caco-2 Permeability Prediction"
task_type: drug_discovery
domain: pharmacokinetics
difficulty: medium
description: biomlbench/tasks/caco2-wang/description.md

data_source:
  type: polaris
  dataset_id: "Caco2_Wang"

dataset:
  answers: caco2-wang/prepared/private/answers.csv
  sample_submission: caco2-wang/prepared/public/sample_submission.csv

grader:
  name: rmse
  grade_fn: biomlbench.tasks.caco2-wang.grade:grade

preparer: biomlbench.tasks.caco2-wang.prepare:prepare

# Biomedical-specific metadata
biomedical_metadata:
  modality: "molecular_properties"
  organ_system: "intestine"
  data_type: "regression"
  clinical_relevance: "drug_absorption"

# Human performance baselines
human_baselines:
  expert_medicinal_chemist: 0.45
  computational_model: 0.52

# Computational requirements
compute_requirements:
  recommended_gpu_memory_gb: 2
  estimated_runtime_minutes: 15
  max_dataset_size_gb: 1
```

## Data Directory Management

The registry supports configurable data directories:

```python
from pathlib import Path
from biomlbench.registry import Registry

# Use custom data directory
custom_registry = Registry(data_dir=Path("/custom/data"))

# Or set on existing registry
new_registry = registry.set_data_dir(Path("/new/data"))

# Get current data directory
data_dir = registry.get_data_dir()
print(f"Data stored in: {data_dir}")
```

## Task Filtering and Organization

### By Biomedical Domain

```python
# Get tasks by specific domains
oncology_tasks = registry.get_tasks_by_domain("oncology")
cardiology_tasks = registry.get_tasks_by_domain("cardiology")
neurology_tasks = registry.get_tasks_by_domain("neurology")

# Available domains
domains = set()
for task_id in registry.list_task_ids():
    task = registry.get_task(task_id)
    domains.add(task.domain)
print(f"Domains: {sorted(domains)}")
```

### By Task Type

```python
# Get tasks by type
imaging_tasks = registry.get_tasks_by_type("medical_imaging")
drug_tasks = registry.get_tasks_by_type("drug_discovery")
protein_tasks = registry.get_tasks_by_type("protein_engineering")

# Available task types
task_types = set()
for task_id in registry.list_task_ids():
    task = registry.get_task(task_id)
    task_types.add(task.task_type)
print(f"Task types: {sorted(task_types)}")
```

### By Difficulty

```python
# Get tasks by difficulty level
easy_tasks = registry.get_tasks_by_difficulty("easy")
medium_tasks = registry.get_tasks_by_difficulty("medium")
hard_tasks = registry.get_tasks_by_difficulty("hard")

# Get lite task set (easy + medium)
lite_tasks = registry.get_lite_task_ids()
```

## Task Metadata Access

### Biomedical Metadata

```python
task = registry.get_task("histopathologic-cancer-detection")

# Access biomedical-specific metadata
metadata = task.biomedical_metadata
print(f"Modality: {metadata['modality']}")
print(f"Organ system: {metadata['organ_system']}")
print(f"Clinical relevance: {metadata['clinical_relevance']}")

# Human baseline performance
if task.human_baselines:
    for role, score in task.human_baselines.items():
        print(f"{role}: {score}")

# Computational requirements
if task.compute_requirements:
    gpu_mem = task.compute_requirements.get('recommended_gpu_memory_gb', 0)
    runtime = task.compute_requirements.get('estimated_runtime_minutes', 0)
    print(f"GPU: {gpu_mem}GB, Runtime: {runtime}min")
```

### Data Source Information

```python
task = registry.get_task("caco2-wang")

# Data source configuration
if task.data_source:
    source_type = task.data_source['type']
    print(f"Data source: {source_type}")
    
    if source_type == "kaggle":
        comp_id = task.data_source['competition_id']
        print(f"Kaggle competition: {comp_id}")
    elif source_type == "polaris":
        dataset_id = task.data_source['dataset_id']
        print(f"Polaris dataset: {dataset_id}")
```

## File Path Resolution

The registry automatically resolves file paths for task components:

```python
task = registry.get_task("caco2-wang")

# Data file paths
print(f"Answers: {task.answers}")
print(f"Sample submission: {task.sample_submission}")

# Directory paths
print(f"Raw data: {task.raw_dir}")
print(f"Public data: {task.public_dir}")
print(f"Private data: {task.private_dir}")

# Verification files
print(f"Checksums: {task.checksums}")
print(f"Leaderboard: {task.leaderboard}")
```

## Task Splits and Collections

### Predefined Splits

```python
# Get tasks from predefined splits
lite_tasks = registry.get_lite_task_ids()

# Check if splits directory exists
splits_dir = registry.get_splits_dir()
if (splits_dir / "medical-imaging.txt").exists():
    with open(splits_dir / "medical-imaging.txt") as f:
        imaging_split = [line.strip() for line in f]
```

### Custom Task Collections

```python
# Create custom task collections
molecular_tasks = []
for task_id in registry.list_task_ids():
    task = registry.get_task(task_id)
    if "molecular" in task.biomedical_metadata.get("modality", ""):
        molecular_tasks.append(task_id)

# Filter by clinical relevance
drug_absorption_tasks = []
for task_id in registry.list_task_ids():
    task = registry.get_task(task_id)
    if task.biomedical_metadata.get("clinical_relevance") == "drug_absorption":
        drug_absorption_tasks.append(task_id)
```

## Error Handling

```python
from biomlbench.registry import registry

try:
    task = registry.get_task("nonexistent-task")
except ValueError as e:
    print(f"Task not found: {e}")

try:
    # Invalid task configuration
    task = registry.get_task("invalid-task")
except ValueError as e:
    print(f"Invalid task config: {e}")
```

## Extending the Registry

### Adding New Tasks

1. Create task directory structure:
   ```bash
   mkdir -p biomlbench/tasks/my-new-task
   ```

2. Add required files:
   - `config.yaml` - Task configuration
   - `prepare.py` - Data preparation logic
   - `grade.py` - Evaluation function
   - `description.md` - Task description

3. Task is automatically discovered on next registry access

### Custom Registry Implementations

```python
from biomlbench.registry import Registry
from pathlib import Path

class CustomRegistry(Registry):
    def get_custom_tasks(self) -> list[str]:
        """Get tasks with custom criteria."""
        return [
            task_id for task_id in self.list_task_ids()
            if self.get_task(task_id).difficulty == "easy"
        ]

# Use custom registry
custom_registry = CustomRegistry(data_dir=Path("/data"))
easy_tasks = custom_registry.get_custom_tasks()
``` 