# Registry API

Task discovery, loading, and organization system for biomedical benchmark tasks.

::: biomlbench.registry
    options:
      show_source: false
      heading_level: 2

## Usage Examples

### Basic Task Access

```python
from biomlbench.registry import registry

# Get task information
task = registry.get_task("caco2-wang")
print(f"Task: {task.name} ({task.domain})")

# List all available tasks
all_tasks = registry.list_task_ids()
print(f"Available tasks: {len(all_tasks)}")
```

### Task Configuration Format

Each task requires a `config.yaml` file:

```yaml
id: caco2-wang
name: "Caco-2 Cell Permeability Prediction"
task_type: drug_discovery
domain: admet
difficulty: medium

data_source:
  type: polaris
  benchmark_id: tdcommons/caco2-wang

dataset:
  answers: caco2-wang/prepared/private/answers.csv
  sample_submission: caco2-wang/prepared/public/sample_submission.csv

grader:
  name: mean-absolute-error
  grade_fn: biomlbench.tasks.caco2-wang.grade:grade

preparer: biomlbench.tasks.caco2-wang.prepare:prepare

biomedical_metadata:
  modality: "molecular"
  data_type: "smiles_regression"
  clinical_relevance: "drug_absorption"

compute_requirements:
  recommended_gpu_memory_gb: 4
  estimated_runtime_minutes: 15
  max_dataset_size_gb: 1
``` 