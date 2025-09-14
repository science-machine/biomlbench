# Installation

This guide covers installing BioML-bench and setting up the environment for running agents on biomedical tasks.

## Prerequisites

- **Python 3.11+**
- **Docker** - For containerized agent execution
- **uv** - Python package manager (recommended)

### Install uv

BioML-bench uses [uv](https://github.com/astral-sh/uv) for fast dependency management:

```bash
# Install uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# Or with pip
pip install uv
```

### Install Docker

Docker is required for secure agent execution:

=== "Ubuntu/Debian"
    ```bash
    sudo apt update
    sudo apt install docker.io
    sudo systemctl start docker
    sudo usermod -aG docker $USER
    # Log out and back in for group changes
    ```

=== "macOS"
    ```bash
    # Install Docker Desktop from https://docker.com/products/docker-desktop
    # Or with Homebrew
    brew install --cask docker
    ```

=== "Windows"
    Download and install [Docker Desktop](https://docker.com/products/docker-desktop).

## Installing BioML-bench

### From Source (Recommended for Development)

```bash
# Clone the repository
git clone https://github.com/science-machine/biomlbench.git
cd biomlbench

# Install dependencies with uv
uv sync

# Activate the virtual environment
source .venv/bin/activate  # Linux/macOS
# or
# .venv\Scripts\activate     # Windows

# Verify installation
biomlbench --help
```

### Development Installation

For contributing to BioML-bench:

```bash
# Clone and install in development mode
git clone https://github.com/science-machine/biomlbench.git
cd biomlbench

# Install with development dependencies
uv sync --extra dev

# Install pre-commit hooks
pre-commit install
```

## Environment Setup

### Option 1: Use Prebuilt Images (Recommended)

The fastest way to get started is to pull our prebuilt Docker images:

```bash
# Pull and tag all prebuilt images
./scripts/pull_prebuilt_images.sh
```

This script pulls the following images and tags them for local use:
- `millerh1/biomlbench-env:v0.1a` → `biomlbench-env`
- `millerh1/aide:v0.1a` → `aide`
- `millerh1/biomni:v0.1a` → `biomni`  
- `millerh1/mlagentbench:v0.1a` → `mlagentbench`
- `millerh1/stella:v0.1a` → `stella`
- `millerh1/dummy:v0.1a` → `dummy`

### Option 2: Build Images Locally

If you prefer to build images locally or need to modify agents:

```bash
# Build the biomlbench-env base image
./scripts/build_base_env.sh

# Build specific agent images
./scripts/build_agent.sh aide
./scripts/build_agent.sh biomni

# Verify the environment
./scripts/test_environment.sh
```


## Configuration

### Polaris API (for polaris-based tasks)

```bash
polaris login --overwrite
```

### Kaggle API (For Kaggle-based Tasks)

Some tasks require Kaggle data. Set up Kaggle API credentials:

```bash
# Download API credentials from https://www.kaggle.com/account
# Place in ~/.kaggle/kaggle.json (Linux/macOS) or %USERPROFILE%\.kaggle\kaggle.json (Windows)

# Set permissions (Linux/macOS only)
chmod 600 ~/.kaggle/kaggle.json
```

### Agent API Keys

For agents that require API access (e.g., AIDE):

```bash
# Create environment file
echo "OPENAI_API_KEY=your-key-here" >> .env
echo "ANTHROPIC_API_KEY=your-key-here" >> .env
echo "GEMINI_API_KEY=your-key-here" >> .env
```

## Verification

Test your installation:

```bash
# Check CLI is working
biomlbench --help

# List available tasks
biomlbench prepare --help

# Test with dummy agent (requires Docker)
biomlbench prepare -t polarishub/tdcommons-caco2-wang
biomlbench run-agent --agent dummy --task-id polarishub/tdcommons-caco2-wang
```


### Getting Help

- Open an issue on [GitHub](https://github.com/science-machine/biomlbench/issues) 