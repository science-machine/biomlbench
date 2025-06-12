# Creating Agents

Guide for developing custom AI agents for BioML-bench.

## Agent Requirements

LLM agents must:

- **Run autonomously** in Docker containers without human intervention
- **Read and comprehend** biomedical task descriptions from `/home/data/description.md`
- **Analyze data** and design appropriate ML approaches
- **Implement and execute** complete ML pipelines (preprocessing, training, evaluation)
- **Generate predictions** and write them to `/home/submission/submission.csv`
- **Handle domain-specific data** (SMILES for molecules, medical images, protein sequences, etc.)

Unlike traditional ML models, agents must solve tasks from scratch using only the task description and training data.

## Agent Structure

```
agents/my-agent/
├── Dockerfile           # Container build instructions
├── config.yaml         # Agent configuration
├── start.sh            # Execution entry point
├── requirements.txt    # Python dependencies
└── src/                # Agent source code
    └── agent.py
```

## Configuration (`config.yaml`)

```yaml
my-agent:
  start: my-agent/start.sh
  dockerfile: my-agent/Dockerfile
  env_vars:
    OPENAI_API_KEY: "${OPENAI_API_KEY}"
  privileged: false  # Set true only if necessary
```

## Docker Setup (`Dockerfile`)

```dockerfile
FROM biomlbench-env

# Copy agent code
COPY src/ /home/agent/src/
COPY requirements.txt /home/agent/

# Install dependencies
WORKDIR /home/agent
RUN /opt/conda/bin/conda run -n agent pip install -r requirements.txt

# Set working directory
WORKDIR /home
```

## Agent Interface

Agents are LLM-based systems that autonomously solve biomedical ML tasks end-to-end. They must:

1. **Understand the task** by reading the description and analyzing the data
2. **Design an ML approach** appropriate for the biomedical problem  
3. **Implement, train, and evaluate** models
4. **Generate predictions** for submission

### Example Agent Structure

```python
# src/agent.py
import pandas as pd
from pathlib import Path

def main():
    # 1. Read and understand the task
    with open('/home/data/description.md', 'r') as f:
        task_description = f.read()
        # Agent parses task type, evaluation metric, clinical context
    
    # 2. Analyze the data
    train_df = pd.read_csv('/home/data/train.csv')
    test_df = pd.read_csv('/home/data/test_features.csv')
    sample_submission = pd.read_csv('/home/data/sample_submission.csv')
    
    # 3. ... Agent Logic ...
    
    # 4. Create submission in required format
    submission = pd.DataFrame({
        'id': test_df['id'],
        'prediction': predictions  # Column name from sample_submission
    })
    
    submission.to_csv('/home/submission/submission.csv', index=False)

if __name__ == "__main__":
    main()
```

### Agent Input Data

- **`description.md`**: Complete task description including:
  - Task type (regression, classification, segmentation)
  - Data format and features
  - Evaluation metric
  - Baseline approaches and references
  
- **`train.csv`**: Training data with features and targets
- **`test_features.csv`**: Test features (no targets)
- **`sample_submission.csv`**: Expected submission format

## Building and Testing

```bash
# Build agent
./scripts/build_agent.sh my-agent # e.g., ./scripts/build_agent.sh aide

# Test agent
biomlbench run-agent --agent my-agent --task-id caco2-wang # e.g., biomlbench run-agent --agent aide --task-id caco2-wang
``` 