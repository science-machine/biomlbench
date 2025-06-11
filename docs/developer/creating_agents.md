# Creating Agents

Guide for developing custom AI agents for BioML-bench.

## Agent Requirements

Agents must:
- Run in Docker containers
- Read task data from `/home/data/`
- Write predictions to `/home/submission/submission.csv`
- Handle biomedical data formats appropriately

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

```python
# src/agent.py
import pandas as pd
from pathlib import Path

def main():
    # Read task description
    with open('/home/data/description.md', 'r') as f:
        description = f.read()
    
    # Load training data
    train_df = pd.read_csv('/home/data/train.csv')
    test_df = pd.read_csv('/home/data/test_features.csv')
    
    # Your ML logic here
    predictions = your_model.predict(test_df)
    
    # Create submission
    submission = pd.DataFrame({
        'id': test_df['id'],
        'prediction': predictions
    })
    
    # Save submission
    submission.to_csv('/home/submission/submission.csv', index=False)

if __name__ == "__main__":
    main()
```

## Building and Testing

```bash
# Build agent
./scripts/build_agent.sh my-agent

# Test agent
biomlbench run-agent --agent my-agent --task-id caco2-wang
``` 