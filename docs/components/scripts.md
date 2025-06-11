# Scripts Reference

The `scripts/` directory contains automation scripts for building, testing, and managing BioML-bench environments and agents.

## Build Scripts

### `build_base_env.sh`

Builds the foundational Docker image (`biomlbench-env`) with biomedical libraries and dependencies.

**Usage:**
```bash
./scripts/build_base_env.sh
```

**Features:**
- Installs Python 3.11 with conda
- Includes biomedical libraries (RDKit, BioPython)
- Configures ML frameworks (TensorFlow, PyTorch)
- Sets up the grading server environment
- Performs post-build validation tests

**Requirements:**
- Docker installed and running
- Internet connection for downloading dependencies
- ~10GB disk space for the image

**Example Output:**
```
🧬 Building BioML-bench Base Environment
=======================================
📋 Pre-build checks...
✅ environment/Dockerfile
✅ environment/requirements.txt
🔨 Building biomlbench-env Docker image...
✅ Successfully built biomlbench-env image
🧪 Testing base image...
✅ Python is available
✅ BioML-bench is importable
✅ Biomedical and ML dependencies are available
✅ Agent conda environment is ready
🎉 Base environment build completed successfully!
```

### `build_agent.sh`

Builds Docker images for individual agents.

**Usage:**
```bash
./scripts/build_agent.sh <agent-name>

# Examples
./scripts/build_agent.sh dummy
./scripts/build_agent.sh aide
```

**Features:**
- Builds agent-specific Docker images
- Inherits from `biomlbench-env` base image
- Configures agent-specific dependencies
- Sets up environment variables and entrypoints

**Agent Directory Structure:**
```
agents/<agent-name>/
├── Dockerfile          # Agent-specific build instructions
├── config.yaml         # Agent configuration
├── start.sh            # Agent execution script
├── requirements.txt    # Additional dependencies
└── src/                # Agent source code
```

### `build_agent_enhanced.sh`

Enhanced agent building with additional features and validation.

**Usage:**
```bash
./scripts/build_agent_enhanced.sh <agent-name> [options]
```

**Additional Features:**
- Multi-stage builds for smaller images
- Build caching optimization
- Security scanning integration
- Performance benchmarking
- Automated testing of built images

## Testing Scripts

### `test_environment.sh`

Comprehensive testing of the BioML-bench environment setup.

**Usage:**
```bash
./scripts/test_environment.sh
```

**Tests Performed:**
- Docker daemon connectivity
- Base image functionality
- Agent image building
- Container execution
- Data mounting and permissions
- Network isolation
- GPU support (if available)

**Example Test Sequence:**
```bash
# Test Docker setup
docker --version
docker info

# Test base environment
docker run --rm biomlbench-env python --version

# Test agent building
./scripts/build_agent.sh dummy

# Test agent execution
biomlbench run-agent --agent dummy --task-id caco2-wang
```

### `validate_environment.sh`

Validates system prerequisites and configuration.

**Usage:**
```bash
./scripts/validate_environment.sh
```

**Validation Checks:**
- Python version compatibility (≥3.11)
- Docker installation and permissions
- Available disk space
- Memory requirements
- GPU support detection
- Network connectivity
- Required system packages

**Prerequisites Validation:**
```bash
# Check Python version
python --version  # Should be ≥3.11

# Check Docker
docker --version
docker run hello-world

# Check disk space
df -h  # Need ~15GB free

# Check memory
free -h  # Recommend ≥8GB RAM
```

### `test_agent_pipeline.py`

End-to-end testing of the agent execution pipeline.

**Usage:**
```bash
python scripts/test_agent_pipeline.py [options]
```

**Test Coverage:**
- Task preparation workflow
- Agent container creation
- Data volume mounting
- Execution environment setup
- Submission file generation
- Cleanup and resource management

**Command Options:**
- `--agent` - Specific agent to test
- `--task` - Specific task to test
- `--timeout` - Execution timeout
- `--retain-containers` - Keep containers for debugging

## Monitoring Scripts

### `monitor_system.py`

System resource monitoring during agent execution.

**Usage:**
```bash
python scripts/monitor_system.py [options]
```

**Monitoring Features:**
- CPU usage tracking
- Memory consumption
- Disk I/O statistics
- GPU utilization (if available)
- Container resource limits
- Network traffic monitoring

**Output Options:**
- Real-time console output
- CSV log files
- JSON metrics export
- Grafana-compatible format

**Example Usage:**
```bash
# Monitor during agent run
python scripts/monitor_system.py --output metrics.csv &
biomlbench run-agent --agent aide --task-id caco2-wang
```

## Utility Scripts

### `robust_postprocess.py`

Post-processing and analysis of agent run results.

**Usage:**
```bash
python scripts/robust_postprocess.py --run-dir <run-directory> [options]
```

**Features:**
- Submission validation and repair
- Log file analysis and summarization
- Performance metric extraction
- Error categorization and reporting
- Resource usage analysis

**Analysis Outputs:**
- Summary reports (JSON/HTML)
- Performance visualizations
- Error logs and diagnostics
- Resource consumption charts

**Command Options:**
- `--run-dir` - Target run directory
- `--output-format` - Report format (json/html/csv)
- `--include-logs` - Include detailed log analysis
- `--validate-submissions` - Check submission format compliance

## Script Configuration

### Environment Variables

Scripts respect these environment variables:

```bash
# Build configuration
export DOCKER_BUILDKIT=1                    # Enable BuildKit
export BUILD_PARALLEL=4                     # Parallel build jobs

# Testing configuration  
export TEST_TIMEOUT=3600                    # Test timeout (seconds)
export RETAIN_TEST_CONTAINERS=false         # Keep test containers

# Resource limits
export MAX_MEMORY="8g"                      # Container memory limit
export MAX_CPUS="4"                         # CPU limit

# Monitoring configuration
export METRICS_INTERVAL=5                   # Monitoring interval (seconds)
export METRICS_OUTPUT="/tmp/metrics"        # Metrics output directory
```

### Build Options

Most build scripts support common options:

```bash
# Build with custom base image
./scripts/build_agent.sh dummy --base-image custom-biomlbench-env

# Skip dependency installation
./scripts/build_agent.sh aide --skip-deps

# Enable debug builds
./scripts/build_agent.sh aide --debug

# Use build cache
./scripts/build_agent.sh aide --cache

# Multi-platform builds
./scripts/build_agent.sh aide --platform linux/amd64,linux/arm64
```

## Common Usage Patterns

### Initial Setup

```bash
# Complete environment setup
./scripts/validate_environment.sh
./scripts/build_base_env.sh
./scripts/test_environment.sh

# Build all agents
for agent in dummy aide mlagentbench; do
    ./scripts/build_agent.sh $agent
done
```

### Development Workflow

```bash
# Build and test specific agent
./scripts/build_agent.sh my-agent --debug
python scripts/test_agent_pipeline.py --agent my-agent --task caco2-wang

# Monitor resource usage
python scripts/monitor_system.py --agent my-agent &
biomlbench run-agent --agent my-agent --task-id caco2-wang

# Analyze results
python scripts/robust_postprocess.py --run-dir runs/latest-run-group/
```

### Continuous Integration

```bash
# CI/CD pipeline
./scripts/validate_environment.sh || exit 1
./scripts/build_base_env.sh || exit 1
./scripts/test_environment.sh || exit 1

# Test all agents
for agent in $(ls agents/*/config.yaml | cut -d/ -f2); do
    ./scripts/build_agent.sh $agent || exit 1
    python scripts/test_agent_pipeline.py --agent $agent --timeout 600 || exit 1
done
```

### Performance Optimization

```bash
# Build with optimizations
export DOCKER_BUILDKIT=1
export BUILD_PARALLEL=8

# Use build cache effectively
./scripts/build_agent.sh aide --cache --parallel

# Monitor and optimize
python scripts/monitor_system.py --profile-memory --profile-cpu
```

## Troubleshooting

### Common Issues

**Build failures:**
```bash
# Clean Docker state
docker system prune -a

# Rebuild base image
./scripts/build_base_env.sh --no-cache

# Check disk space
df -h
```

**Test failures:**
```bash
# Run validation first
./scripts/validate_environment.sh

# Check Docker permissions
docker run hello-world

# Verify Python environment
python -c "import biomlbench; print('OK')"
```

**Resource issues:**
```bash
# Monitor system resources
python scripts/monitor_system.py --real-time

# Adjust container limits
export MAX_MEMORY="16g"
export MAX_CPUS="8"
```

### Script Dependencies

Most scripts require:
- Bash shell (Linux/macOS)
- Docker (version ≥20.10)
- Python ≥3.11
- Internet connectivity
- Sufficient disk space (≥15GB)

For GPU support:
- NVIDIA drivers
- NVIDIA Container Toolkit
- CUDA-compatible hardware 