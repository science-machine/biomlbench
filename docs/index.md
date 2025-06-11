# BioML-bench

A benchmark suite for evaluating machine learning agents on biomedical tasks.

BioML-bench is built on top of [MLE-bench](https://github.com/openai/mle-bench) and provides a comprehensive framework for benchmarking AI agents on biomedical machine learning (BioML) tasks including protein engineering, drug discovery, medical imaging, and clinical biomarkers.

## 🧬 Features

- **Diverse Biomedical Tasks**: Protein engineering, drug discovery, medical imaging, clinical biomarkers
- **Agent-Agnostic Evaluation**: Any agent that can produce CSV outputs can be evaluated
- **Standardized Metrics**: Domain-specific evaluation metrics (RMSD, TM-score, AUC-ROC, etc.)
- **Human Baselines**: Built-in human performance benchmarks for comparison
- **Secure Evaluation**: Containerized execution with no data leakage
- **Extensible Framework**: Easy to add new biomedical tasks
- **Biomedical Libraries**: Pre-installed RDKit, BioPython, and other domain-specific tools

## 🚀 Quick Start

```bash
# Install with uv
git clone https://github.com/bioml-bench/bioml-bench.git
cd bioml-bench
uv sync

# Prepare a task
biomlbench prepare -t caco2-wang

# Run an agent
biomlbench run-agent --agent dummy --task-id caco2-wang

# Grade results (submission.jsonl is auto-generated)
biomlbench grade --submission <run-group-dir>/submission.jsonl --output-dir results/
```

## 📊 Available Tasks

### Medical Imaging
- **histopathologic-cancer-detection**: Cancer detection in histopathology patches

### Drug Discovery  
- **caco2-wang**: Molecular property prediction (intestinal permeability)

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