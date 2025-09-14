# BioML-bench

A benchmark suite for evaluating LLM agents on biomedical machine learning tasks.

![BioML-bench Overview](images/biombench_ga.png)

**📄 Paper**: [BioML-bench: Evaluation of AI Agents for End-to-End Biomedical ML](https://www.biorxiv.org/content/10.1101/2025.09.01.673319v2)

BioML-bench is built on top of [MLE-bench](https://github.com/openai/mle-bench) and provides a comprehensive framework for benchmarking LLM agents on biomedical machine learning (BioML) tasks including protein engineering, drug discovery, single cell omics, medical imaging, and clinical biomarkers.

Agents autonomously read task descriptions, analyze biomedical data, design appropriate ML approaches, and implement complete solutions from scratch.

## 🧬 Features

- **Diverse Biomedical Tasks**: Protein engineering, drug discovery, single cell omics, medical imaging, clinical biomarkers
- **Agent-Agnostic Evaluation**: Any LLM agent that can read task descriptions and produce file/folder submissions can be evaluated
- **Human Baselines**: Built-in human performance benchmarks for comparison
- **Extensible Framework**: Easy to add new biomedical tasks
- **Biomedical Libraries**: Pre-installed RDKit, BioPython, and other domain-specific tools for use by agents

## 🚀 Quick Start

Install the package (requires [uv](https://docs.astral.sh/uv/)):

```bash
# Install with uv
git clone https://github.com/science-machine/biomlbench.git
cd biomlbench
uv sync

# Activate the environment
source .venv/bin/activate  # Linux/macOS
# or .venv\Scripts\activate  # Windows
```

Pull prebuilt Docker images and run a benchmark:

```bash
# Pull prebuilt images (recommended - saves build time)
./scripts/pull_prebuilt_images.sh

# Prepare a task
biomlbench prepare -t polarishub/tdcommons-caco2-wang

# Run an agent (in this case, a dummy agent that returns null predictions)
biomlbench run-agent --agent dummy --task-id polarishub/tdcommons-caco2-wang

# Grade results (submission.jsonl is auto-generated)
biomlbench grade --submission <run-group-dir>/submission.jsonl --output-dir results/
```

### Full v0.1a Benchmark

To evaluate on our complete 24-task benchmark from the preprint:

```bash
# Make sure to have a .env file with the following variables:
echo "OPENAI_API_KEY=<your_openai_api_key>" >> .env
echo "ANTHROPIC_API_KEY=<your_anthropic_api_key>" >> .env
echo "OPENROUTER_API_KEY=<your_openrouter_api_key>" >> .env
echo "GEMINI_API_KEY=<your_gemini_api_key>" >> .env
echo "MEM0_API_KEY=<your_mem0_api_key>" >> .env

# Pull prebuilt images
./scripts/pull_prebuilt_images.sh

# Prepare all v0.1a tasks (downloads several GB of data)
./scripts/prepare_tasks_from_file.sh experiments/biomlbench_v0.1a.txt

# Run all agents on the full benchmark
biomlbench run-agent --agent dummy --task-list experiments/biomlbench_v0.1a.txt
biomlbench run-agent --agent aide --task-list experiments/biomlbench_v0.1a.txt
biomlbench run-agent --agent biomni --task-list experiments/biomlbench_v0.1a.txt
biomlbench run-agent --agent mlagentbench --task-list experiments/biomlbench_v0.1a.txt
biomlbench run-agent --agent stella --task-list experiments/biomlbench_v0.1a.txt
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
- **Tasks** (`biomlbench/tasks/`) - Individual biomedical benchmark tasks, each of which contains at least 1 dataset
- **Scripts** (`scripts/`) - Build, test, and deployment automation
- **Deploy** (`deploy/`) - Deployment scripts for running BioML-bench on GCP and AWS (Under development)

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