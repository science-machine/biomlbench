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
