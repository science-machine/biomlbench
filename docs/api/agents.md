# Agents API

The agents module provides functionality for running AI agents in containerized environments on biomedical tasks.

## Core Functions

::: biomlbench.agents.run_agent_async

::: biomlbench.agents.run_agent

## Data Classes

::: biomlbench.agents.AgentTask

## Worker Functions

::: biomlbench.agents.worker

## Utility Functions

::: biomlbench.agents.get_task_ids

::: biomlbench.agents.create_task_list_file

## Agent Execution Workflow

```python
import asyncio
from biomlbench.agents import run_agent_async

# Run agent asynchronously
run_group, submission_path = await run_agent_async(
    agent_id="dummy",
    task_ids=["caco2-wang"],
    n_workers=1,
    n_seeds=1
)

print(f"Run completed: {run_group}")
print(f"Submission ready: {submission_path}")
```

## Container Management

The agents module handles Docker container lifecycle:

1. **Container Creation** - Spins up task-specific containers
2. **Volume Mounting** - Mounts task data and output directories
3. **Environment Setup** - Configures environment variables
4. **Agent Execution** - Runs agent start script
5. **Output Collection** - Extracts submissions and logs
6. **Cleanup** - Removes containers (unless retained)

## Parallel Execution

Support for running multiple tasks in parallel:

```python
# Run multiple tasks with parallel workers
await run_agent_async(
    agent_id="aide",
    task_ids=["caco2-wang", "histopathologic-cancer-detection"],
    n_workers=4,  # 4 parallel workers
    n_seeds=2     # 2 random seeds per task
)
```

## Configuration

Agent execution can be customized via:

- **Container config** - Docker resource limits and settings
- **Environment variables** - Agent-specific configuration
- **Data directory** - Custom data storage location
- **Retention policy** - Keep containers for debugging 