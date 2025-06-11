# Architecture Overview

BioML-bench follows a modular architecture designed for extensibility and security.

## Core Components

### Registry System
- **Tasks** - Discovery and management of biomedical benchmark tasks
- **Agents** - Registration and configuration of AI agents
- **Data Sources** - Pluggable data acquisition from multiple repositories

### Execution Engine
- **Containerization** - Docker-based secure execution environment
- **Parallel Processing** - Multi-worker task execution
- **Resource Management** - Memory, CPU, and GPU allocation

### Evaluation Framework
- **Domain-Specific Metrics** - Biomedical evaluation functions
- **Human Baselines** - Expert performance comparisons
- **Medal System** - Kaggle-style performance ranking

## Security Model

### Container Isolation
- **Network isolation** - No external network access during execution
- **Filesystem isolation** - Read-only data access, private answer protection
- **User isolation** - Non-root execution, limited capabilities

### Data Protection
- **Private/public splits** - Clear separation of training and evaluation data
- **Checksum validation** - Data integrity verification
- **Access controls** - Restricted access to ground truth labels

## Data Flow

1. **Task Preparation** - Download and process datasets
2. **Agent Execution** - Run agents in isolated containers
3. **Submission Collection** - Extract predictions and logs
4. **Evaluation** - Score submissions against private answers
5. **Reporting** - Generate comprehensive evaluation reports

## Extension Points

### Adding New Components
- **Tasks** - Implement prepare/grade functions
- **Agents** - Create Docker containers with standard interface
- **Data Sources** - Implement DataSource interface
- **Metrics** - Define custom evaluation functions

### Configuration System
- **YAML-based** - Human-readable configuration files
- **Environment variables** - Runtime configuration options
- **Container configs** - Docker resource and security settings 