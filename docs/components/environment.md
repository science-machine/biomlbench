# Environment Reference

The `environment/` directory contains the containerized execution environment for BioML-bench agents, including the base Docker image, grading server, and configuration files.

## Overview

BioML-bench uses Docker containers to provide secure, isolated execution environments for AI agents. The environment includes:

- **Base Docker image** with biomedical libraries
- **Grading server** for submission validation
- **Container configuration** system
- **Security isolation** mechanisms

## Core Components

### Base Docker Image (`Dockerfile`)

The foundational Docker image provides a complete biomedical ML environment:

**Base System:**
- Ubuntu 22.04 LTS
- Python 3.11 via Miniconda
- Essential system packages and development tools

**Biomedical Libraries:**
- **RDKit** - Molecular informatics toolkit
- **BioPython** - Biological computation library
- **scikit-learn** - Machine learning framework
- **pandas, numpy** - Data processing

**ML Frameworks:**
- **TensorFlow 2.17** with CUDA support
- **PyTorch 2.2.0** with GPU acceleration
- **Transformers** for NLP models
- **OpenCV** for computer vision

**Environment Structure:**
```dockerfile
FROM ubuntu:22.04

# System packages and Python
RUN apt-get update && apt-get install -y \
    python3 python3-pip python3-venv \
    build-essential git curl wget \
    libsm6 libxext6 ffmpeg

# Miniconda installation
RUN wget miniconda.sh && bash miniconda.sh -b -p /opt/conda

# Create isolated environments
RUN conda create -n agent python=3.11     # Agent execution
RUN conda create -n biomlb python=3.11    # Grading server

# Install biomedical dependencies
RUN conda install -n agent -c conda-forge \
    rdkit biopython scikit-learn pandas numpy
```

### Grading Server (`grading_server.py`)

A Flask-based validation server that runs inside containers to evaluate agent submissions:

::: environment.grading_server

**Key Features:**
- HTTP health check endpoint (`/health`)
- Submission validation against private answers
- Real-time evaluation feedback
- Resource monitoring and limits

**Usage in Containers:**
```python
# Server starts automatically via entrypoint.sh
# Agents can check health status
curl http://localhost:5000/health

# Server validates submissions internally
# No direct agent interaction required
```

### Container Entrypoint (`entrypoint.sh`)

The entry script that configures the container environment and starts services:

**Responsibilities:**
- User environment setup (non-root execution)
- Directory permissions configuration
- Grading server initialization
- Agent execution orchestration

**Execution Flow:**
1. Validate container environment
2. Set up user permissions and directories
3. Start grading server in background
4. Wait for server health check
5. Execute agent start script
6. Clean up resources on exit

**Security Features:**
- Non-root user execution (`nonroot`)
- Private directory isolation (`/private/`)
- Read-only data mounts
- Resource limit enforcement

### Task Instructions

Standardized instructions provided to agents:

#### `instructions.txt`
Complete task instructions with full context and examples.

#### `instructions_obfuscated.txt`  
Minimal instructions to prevent overfitting and data leakage.

**Instruction Format:**
```markdown
# Task: [Task Name]

## Objective
[Clear description of the biomedical problem]

## Data Description
- train.csv: Training data with features and targets
- test_features.csv: Test data (features only)
- sample_submission.csv: Expected submission format

## Evaluation Metric
[Domain-specific metric: AUC-ROC, RMSE, etc.]

## Submission Format
CSV file with columns: [id, prediction]

## Important Notes
- No external data sources allowed
- Use only provided training data
- Submit predictions for all test samples
```

### Validation Script (`validate_submission.sh`)

Shell script for basic submission format validation:

```bash
#!/bin/bash
# Validate submission file exists and has correct format

SUBMISSION_FILE="/home/submission/submission.csv"

if [[ ! -f "$SUBMISSION_FILE" ]]; then
    echo "ERROR: submission.csv not found"
    exit 1
fi

# Check CSV format
if ! head -1 "$SUBMISSION_FILE" | grep -q "id.*prediction"; then
    echo "ERROR: Invalid CSV headers"
    exit 1
fi

echo "Submission validation passed"
```

## Container Configuration

### Default Configuration (`config/container_configs/default.json`)

```json
{
    "mem_limit": null,        # No memory limit (use system default)
    "shm_size": "4G",        # Shared memory for large datasets
    "nano_cpus": 4000000000  # 4 CPU cores
}
```

### Custom Configuration Options

**Resource Limits:**
```json
{
    "mem_limit": "8g",           # 8GB memory limit
    "shm_size": "4g",           # 4GB shared memory
    "nano_cpus": 8000000000,    # 8 CPU cores
    "gpus": -1,                 # All available GPUs
    "runtime": "sysbox-runc"    # Enhanced security runtime
}
```

**Security Settings:**
```json
{
    "privileged": false,         # Disable privileged mode
    "user": "nonroot",          # Non-root user execution
    "read_only": false,         # Allow container writes
    "security_opt": [           # Security options
        "no-new-privileges:true"
    ]
}
```

**Volume Mounting:**
```json
{
    "volumes": {
        "/host/data": {
            "bind": "/home/data",
            "mode": "ro"            # Read-only data access
        },
        "/host/output": {
            "bind": "/home/submission",
            "mode": "rw"            # Write access for submissions
        }
    }
}
```

## Environment Variables

### System Environment

Set automatically by the container:

```bash
# Directory paths
DATA_DIR="/home/data"
SUBMISSION_DIR="/home/submission"
LOGS_DIR="/home/logs"
CODE_DIR="/home/code"
AGENT_DIR="/home/agent"

# Python environment
PYTHONPATH="/opt/conda/envs/agent/lib/python3.11/site-packages"
CONDA_DEFAULT_ENV="agent"

# Task information
TASK_ID="caco2-wang"            # Set by agent runner
TASK_TYPE="drug_discovery"      # Biomedical task type
```

### Agent-Specific Environment

Configurable via agent configuration:

```bash
# API keys (for agents that need them)
OPENAI_API_KEY="sk-..."
ANTHROPIC_API_KEY="sk-ant-..."

# Agent configuration
AGENT_TIMEOUT="3600"            # 1 hour timeout
MAX_ITERATIONS="10"             # Iteration limit
VERBOSE="true"                  # Enable verbose logging
```

## Security Model

### Isolation Mechanisms

**Container Isolation:**
- Separate network namespace
- Isolated filesystem
- Resource limits and quotas
- User namespace isolation

**Data Protection:**
- Private answers in protected directory (`/private/`)
- Read-only access to task data
- No network access during execution
- Temporary filesystem for scratch space

**Process Security:**
- Non-root user execution
- Limited capabilities
- No privileged operations
- Process monitoring and limits

### Sysbox Runtime (Recommended)

Enhanced container security using [Sysbox](https://github.com/nestybox/sysbox):

**Benefits:**
- VM-like isolation
- Secure nested containers
- Better resource management
- Reduced attack surface

**Configuration:**
```json
{
    "runtime": "sysbox-runc",
    "security_opt": [
        "systempaths=unconfined"
    ]
}
```

### Privileged Containers (Special Cases)

Some agents (e.g., OpenHands) require privileged access:

**Requirements:**
- Set `I_ACCEPT_RUNNING_PRIVILEGED_CONTAINERS=True`
- Understand security implications
- Use only in trusted environments

**Configuration:**
```yaml
# agents/opendevin/config.yaml
opendevin:
  privileged: true
  start: opendevin/start.sh
  dockerfile: opendevin/Dockerfile
```

## Volume Mounting Strategy

### Data Volume Structure

```
/home/data/                     # Task data (read-only)
├── description.md              # Task description
├── train.csv                  # Training data
├── test_features.csv          # Test features
├── sample_submission.csv      # Expected format
└── human_baselines.csv        # Human performance (if available)

/home/submission/              # Agent output (read-write)
└── submission.csv             # Agent predictions

/private/data/task-id/         # Private evaluation data
└── prepared/private/
    └── answers.csv            # Ground truth (inaccessible to agents)
```

### Mount Configuration

Configured by the agent runner:

```python
volumes_config = {
    # Task data (read-only)
    str(task.public_dir): {
        "bind": "/home/data",
        "mode": "ro"
    },
    # Private answers (read-only, restricted access)
    str(task.private_dir): {
        "bind": f"/private/data/{task.id}/prepared/private/",
        "mode": "ro"
    }
}
```

## Build Process

### Multi-Stage Build

The Dockerfile uses multi-stage builds for efficiency:

```dockerfile
# Stage 1: Base system and conda
FROM ubuntu:22.04 as base
RUN apt-get update && install system packages
RUN install miniconda

# Stage 2: Biomedical dependencies  
FROM base as biomedical
RUN conda install biomedical libraries
RUN pip install ML frameworks

# Stage 3: Final image
FROM biomedical as final
COPY application files
RUN final configuration
```

### Build Optimization

**Caching Strategy:**
- Layer-wise dependency installation
- Separate system and Python packages
- Conditional heavy dependency installation

**Size Optimization:**
- Multi-stage builds
- Cleanup of package managers
- Removal of development headers

**Example Build Command:**
```bash
docker build \
    --target final \
    --build-arg INSTALL_HEAVY_DEPENDENCIES=true \
    --tag biomlbench-env \
    -f environment/Dockerfile \
    .
```

## Troubleshooting

### Common Issues

**Container fails to start:**
```bash
# Check Docker daemon
docker info

# Verify image exists
docker images biomlbench-env

# Check container logs
docker logs <container-id>
```

**Permission denied errors:**
```bash
# Check volume permissions
ls -la /path/to/mounted/directory

# Verify user mapping
docker exec <container> id

# Check SELinux context (Linux)
ls -Z /path/to/mounted/directory
```

**Grading server not responding:**
```bash
# Check server health
docker exec <container> curl http://localhost:5000/health

# Verify server process
docker exec <container> ps aux | grep grading_server

# Check server logs
docker exec <container> cat /var/log/grading_server.log
```

### Resource Monitoring

**Memory usage:**
```bash
# Monitor container memory
docker stats <container-id>

# Check OOM kills
dmesg | grep -i "killed process"
```

**Disk usage:**
```bash
# Container filesystem usage
docker exec <container> df -h

# Docker space usage
docker system df
```

**GPU usage (if applicable):**
```bash
# Check GPU availability
docker exec <container> nvidia-smi

# Monitor GPU utilization
nvidia-smi -l 1
``` 