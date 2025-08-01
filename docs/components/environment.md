# Environment Reference

The `environment/` directory contains the containerized execution environment for BioML-bench agents, including the base Docker image, grading server, and configuration files.

## Overview

BioML-bench uses Docker containers to provide secure, isolated execution environments for AI agents. The environment includes:

- **Base Docker image** with biomedical libraries
- **Grading server** for submission validation
- **Container configuration** system

## Core Components

### Base Docker Image (`Dockerfile`)

The foundational Docker image provides a complete biomedical ML environment:

**Base System:**

- Ubuntu 22.04 LTS
- Python 3.11 via Miniconda
- Essential system packages and development tools

**Libraries:**

- **RDKit** - Molecular informatics toolkit
- **BioPython** - Biological computation library
- **scikit-learn** - Machine learning framework
- **pandas, numpy** - Data processing
- **TensorFlow** with CUDA support
- **PyTorch** with GPU acceleration


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

Standardized instructions provided to agents. Includes info about general task structure, data, validation, and submission details.

`instructions.txt` - Complete task instructions.

`instructions_obfuscated.txt`  (TODO: NEEDS TO BE CHECKED) - Minimal instructions to prevent overfitting and data leakage.

### Validation Script (`validate_submission.sh`)

Shell script for basic submission format validation. Agents are instructed to run this script to validate their submissions before finishing.

## Container Configuration

### Default Configuration (`config/container_configs/default.json`)

```json
{
    "mem_limit": null,       # No memory limit (use system default)
    "shm_size": "4G",        # Shared memory for large datasets
    "nano_cpus": 4e9         # 4 CPU cores
}
```

### Custom Configuration Options

**Resource Limits:**
```json
{
    "mem_limit": "8g",          # 8GB memory limit
    "shm_size": "4g",           # 4GB shared memory
    "nano_cpus": 8e9,           # 8 CPU cores
    "gpus": -1,                 # All available GPUs
    "runtime": "sysbox-runc"    # Enhanced security runtime
}
```

**Security Settings:**
```json
{
    "privileged": false,        # Disable privileged mode
    "user": "nonroot",          # Non-root user execution
    "read_only": false,         # Allow container writes
    "security_opt": [           # Security options
        "no-new-privileges:true"
    ]
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
```


## Volume Mounting Strategy

### Data Volume Structure

```
/home/data/                    # Task data (read-only)
├── description.md             # Task description
├── train.<ext>                # Training data (e.g., train.csv, train.h5ad)
├── test_features.<ext>        # Test features (e.g., test_features.csv, test_features.h5ad)
├── sample_submission.<ext>    # Expected format (e.g., sample_submission.csv, sample_submission.h5ad)
└── human_baselines.<ext>      # Human performance (if available) (e.g., human_baselines.csv, human_baselines.h5ad)

/home/submission/              # Agent output (read-write)
└── submission.<ext>           # Agent predictions (e.g., submission.csv, submission.h5ad)

/private/data/task-id/         # Private evaluation data
└── prepared/private/
    └── answers.<ext>          # Ground truth (inaccessible to agents) (e.g., answers.csv, answers.h5ad)
```


## Build Process



**Size Optimization:**
- Multi-stage builds
- Cleanup of package managers
- Removal of development headers

**Example Build Command (NEEDS TO BE CHECKED):**
```bash
docker build \
    --target final \
    --build-arg INSTALL_HEAVY_DEPENDENCIES=true \
    --tag biomlbench-env \
    -f environment/Dockerfile \
    .
```
