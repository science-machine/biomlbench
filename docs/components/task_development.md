# Task Development

Guide for creating new biomedical benchmark tasks in BioML-bench.

## Task Structure

Each task requires:

- `config.yaml` - Task configuration and metadata
- `prepare.py` - Data preparation logic
- `grade.py` - Evaluation function
- `description.md` - Task description for agents

## Quick Start

1. Create task directory:
   ```bash
   mkdir -p biomlbench/tasks/my-new-task
   ```

2. Add configuration file
3. Implement preparation logic
4. Define evaluation metrics
5. Test with dummy agent

See [Adding Tasks Guide](../developer/adding_tasks.md) for detailed instructions. 