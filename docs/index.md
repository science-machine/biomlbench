# BioML-bench

A benchmark suite for evaluating LLM agents on biomedical machine learning tasks.

BioML-bench is built on top of [MLE-bench](https://github.com/openai/mle-bench) and provides a comprehensive framework for benchmarking LLM agents on biomedical machine learning (BioML) tasks including protein engineering, drug discovery, medical imaging, and clinical biomarkers.

Agents autonomously read task descriptions, analyze biomedical data, design appropriate ML approaches, and implement complete solutions from scratch.

## 🧬 Features

- **Diverse Biomedical Tasks**: Protein engineering, drug discovery, medical imaging, clinical biomarkers
- **Agent-Agnostic Evaluation**: Any LLM agent that can read task descriptions and produce CSV submissions can be evaluated
- **Human Baselines**: Built-in human performance benchmarks for comparison
- **Secure Evaluation**: Containerized execution with no data leakage
- **Extensible Framework**: Easy to add new biomedical tasks
- **Biomedical Libraries**: Pre-installed RDKit, BioPython, and other domain-specific tools

## 🚀 Quick Start

Install the package (requires [uv](https://docs.astral.sh/uv/)):

```bash
# Install with uv
git clone https://github.com/science-machine/biomlbench.git
cd biomlbench
uv sync
```

Build and run a benchmark with an agent:

```bash
# Prepare a task
biomlbench prepare -t polarishub/tdcommons-caco2-wang

# Run an agent (in this case, a dummy agent that returns null predictions)
biomlbench run-agent --agent dummy --task-id polarishub/tdcommons-caco2-wang

# Grade results (submission.jsonl is auto-generated)
biomlbench grade --submission <run-group-dir>/submission.jsonl --output-dir results/
```

## 📊 Available Tasks

### Medical Imaging
- **manual/histopathologic-cancer-detection**: Cancer detection in histopathology patches

### Drug Discovery  
- **polarishub/tdcommons-caco2-wang**: Molecular property prediction (intestinal permeability)
- **polarishub/*** : 80+ drug discovery and molecular property prediction tasks from Polaris

## 🏗️ Architecture

BioML-bench follows a modular architecture:

- **Core Framework** (`biomlbench/`) - Task management, grading, data handling
- **Agent System** (`agents/`) - Agent registry and execution framework  
- **Environment** (`environment/`) - Containerized execution environment
- **Tasks** (`biomlbench/tasks/`) - Individual biomedical benchmark tasks
- **Scripts** (`scripts/`) - Build, test, and deployment automation

## 📚 Documentation Structure

- **[Quick Start](installation.md)** - Get up and running
- **[API Reference](api/overview.md)** - Complete API documentation
- **[Components](components/environment.md)** - Deep dive into system components
- **[Developer Guide](developer/contributing.md)** - Extend the framework

## 🤝 Contributing

We welcome contributions! See our [Contributing Guide](developer/contributing.md) for details on:

- Adding new biomedical tasks
- Creating custom agents
- Extending data sources
- Improving documentation 