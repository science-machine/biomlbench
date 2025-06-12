# BioML-bench

**THIS WAS GENERATED IN COLLABORATION WITH CLAUDE. EXPECT BUGS.**

**A benchmark suite for evaluating machine learning agents on biomedical tasks**

BioML-bench is built on top of [MLE-bench](https://github.com/openai/mle-bench) and provides a comprehensive framework for benchmarking AI agents on biomedical machine learning (BioML) tasks including protein engineering, drug discovery, medical imaging, and clinical biomarkers.

## Development TODOs (unordered)

- [ ] Add `experiments/splits` which will contain files like `medical-imaging.txt` that contain a list of task ids for only a given domain, or type of benchmarking run (e.g., lite, all).
- [ ] Decide if we want to have `experiments/familiarity` like the original MLE-bench. This would involve measuring the relationship between model performance and the "familiarity" of the model for the task description. This tries to assess if the model is "cheating" by already understanding the task.
- [ ] Need to decide if we can add human baselines in for Polaris. We would probably have to manually scrape the leaderboard (possibly with playwright) and then we would need to add a `supports_human_baselines` method to the `PolarisDataSource` class which returns `True`. And we would need to add a `get_human_baselines` method to the `PolarisDataSource` class which returns a pandas dataframe with the human baselines.
- [ ] Need an automated pipeline for scraping all the polaris and kaggle tasks and creating tasks/ directories for them.
- [ ] Revisit the `generate_submission` method. It currently assumes that the targets in train.csv are the last column. It also requires data in CSV format...
- [ ] Update the models used by the AIDE agent.
- [ ] Find a way to provide task-specific Docker containers for agents. Additionally, need to figure out whether we can enable more flexible base environment container (not hardcoding pytorch version, etc.)
- [ ] Figure out what they did in [the MLE-bench version of aideml](https://github.com/WecoAI/aideml/compare/main...thesofakillers:aideml:main) and see if we should implemnt that.
- [ ] Need to figure out what to do when no human baselines are available. How would scoring work? Also can we realistically combine scores across databases of benchmarks?
- [ ] Need to decide whether we want to incude obfuscated instructions/descriptions for tasks. This would probably take a lot of work but would be useful for benchmarking.


## How to wrap a new benchmark database

Example: Polaris.

To wrap polaris, we needed to:

1. Understand the data source. In this case, Polaris has a high-level API that can be used to automatically download and split data.
2. Create a new data source module in `biomlbench/data_sources/`. In the case of Polaris, we created `biomlbench/data_sources/polaris.py`.
3. Implement the `XYZDataSource` class. This class must implement the `download` and `get_leaderboard` methods.

## 🧬 Features

- **Diverse Biomedical Tasks**: Protein engineering, drug discovery, medical imaging, clinical biomarkers
- **Agent-Agnostic Evaluation**: Any agent that can produce CSV outputs can be evaluated
- **Human Baselines**: Built-in human performance benchmarks for comparison
- **Secure Evaluation**: Containerized execution with no data leakage
- **Extensible Framework**: Easy to add new biomedical tasks
- **Biomedical Libraries**: Pre-installed RDKit, BioPython, and other domain-specific tools

## 🚀 Quick Start

### Installation

BioML-bench uses [uv](https://github.com/astral-sh/uv) for dependency management:

```bash
# Install uv if you haven't already
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone the repository
git clone https://github.com/science-machine/biomlbench.git
cd biomlbench

# Install dependencies
uv sync
```

### Environment Setup

Before running agents, build the base Docker environment:

```bash
# Build the base environment with biomedical libraries
./scripts/build_base_env.sh

# Verify the environment is working
./scripts/test_environment.sh
```

### Basic Usage

1. **Prepare a task dataset:**
```bash
biomlbench prepare -t caco2-wang
```

2. **Run baselines to generate submissions:**
```bash
biomlbench run-baseline --baseline all caco2-wang
biomlbench grade --submission submission.jsonl --output-dir results/
```

3. **Run an agent on a task:**

```bash
biomlbench run-agent --agent dummy --task-id caco2-wang 
# The submission file is automatically generated in the run directory
biomlbench grade --submission <run-group-directory>/submission.jsonl --output-dir results/
```

4. **Grade a specific submission (optional)**

This is useful for manual/external submissions.

```bash
biomlbench grade-sample /path/to/submission.csv caco2-wang
```

## 📊 Available Tasks

### Medical Imaging
- `histopathologic-cancer-detection`: Cancer detection in histopathology patches
- *More tasks coming soon...*

### Protein Engineering
- *Coming soon: ProteinGym fitness prediction, structure prediction*

### Drug Discovery  
- `caco2-wang`: Molecular property prediction (intestinal permeability)
- *More tasks coming soon...*

### Clinical Biomarkers
- *Coming soon: Biomarker discovery, clinical outcome prediction*

## 🏗️ Task Structure

Each task in BioML-bench follows a standardized structure:

```
biomlbench/tasks/task-name/
├── config.yaml              # Task metadata and configuration
├── prepare.py               # Data preparation script
├── grade.py                 # Evaluation logic
├── description.md           # Full task description for agents
├── description_obfuscated.md # Minimal description (anti-overfitting)
├── checksums.yaml           # Data integrity verification
└── leaderboard.csv          # Historical performance baselines
```

### Task Configuration

Tasks are configured via `config.yaml`:

```yaml
id: task-name
name: "Human-readable Task Name"
task_type: medical_imaging  # or protein_engineering, drug_discovery, etc.
domain: oncology            # specific biomedical domain
difficulty: medium          # easy, medium, hard
description: path/to/description.md

dataset:
  answers: path/to/private/answers.csv
  sample_submission: path/to/public/sample_submission.csv

grader:
  name: auc-roc
  grade_fn: biomlbench.tasks.task-name.grade:grade

preparer: biomlbench.tasks.task-name.prepare:prepare

# Biomedical-specific metadata
biomedical_metadata:
  modality: "histopathology"
  organ_system: "breast"
  data_type: "image_classification"
  clinical_relevance: "cancer_detection"

# Human baselines (when available)
human_baselines:
  expert_physician: 0.95
  resident_physician: 0.87

# Computational requirements
compute_requirements:
  recommended_gpu_memory_gb: 8
  estimated_runtime_minutes: 30
  max_dataset_size_gb: 8
```

## 🤖 Agent Evaluation

### Environment Requirements

Agents run in Docker containers based on the `biomlbench-env` image, which includes:
- Python 3.11 with scientific computing libraries
- RDKit for molecular data processing
- BioPython for biological sequence analysis
- Common ML frameworks (scikit-learn, pandas, numpy)
- Biomedical-specific tools and dependencies

### For Agent Developers

To evaluate your agent on BioML-bench:

#### 1. Create Agent Docker Image

Create a `Dockerfile` for your agent:

```dockerfile
FROM biomlbench-env

# Copy your agent code
COPY agent.py /home/agent/
COPY requirements.txt /home/agent/

# Install any additional dependencies
WORKDIR /home/agent
RUN pip install -r requirements.txt

# Set the entrypoint
CMD ["python", "/home/agent/agent.py"]
```

Build your agent:
```bash
docker build -t my-agent .
```

#### 2. Agent Interface

Your agent should follow this interface:

**Input locations:**
- Task description: `/home/data/description.md`
- Training data: `/home/data/train.csv` 
- Test features: `/home/data/test_features.csv`
- Sample submission: `/home/data/sample_submission.csv`

**Output location:**
- Predictions: `/home/submission/submission.csv`

**Agent code structure:**
```python
import pandas as pd
import os

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
        'id': test_df['id'],  # or appropriate ID column
        'prediction': predictions
    })
    
    # Save submission
    submission.to_csv('/home/submission/submission.csv', index=False)

if __name__ == "__main__":
    main()
```

#### 3. Run Agent Evaluation

```bash
# Run agent on a specific task
biomlbench run-agent --agent my-agent --task-id caco2-wang

# Run agent on multiple tasks
biomlbench run-agent --agent my-agent --task-id histopathologic-cancer-detection
```

#### 4. Process Results

```bash
# Submission file is automatically generated during agent execution
# Grade the submission directly
biomlbench grade --submission <run-group-directory>/submission.jsonl --output-dir results/
```

### Built-in Agents

BioML-bench includes several reference agents:

- **Dummy Agent**: Simple baseline that generates random predictions
- **AIDE Agent**: Automated research agent that can solve ML tasks end-to-end

To use AIDE agent (requires OpenAI API key):
```bash
# Set up API key
echo "OPENAI_API_KEY=your-key-here" > .env

# Run AIDE agent
biomlbench run-agent --agent aide --task-id caco2-wang
```

### Agent Execution Environment

- **Security**: Agents run in isolated containers with no network access during execution
- **Resources**: Configurable CPU/memory limits per task requirements  
- **Timeout**: Agents have task-specific time limits to prevent infinite runs
- **Biomedical Libraries**: Pre-installed domain-specific tools (RDKit for molecules, BioPython for sequences)

### Task Splits

Use predefined task splits for standardized evaluation:

```bash
# Medical imaging tasks only
biomlbench run-agent --agent my-agent --task-list experiments/splits/medical-imaging.txt

# Drug discovery tasks only  
biomlbench run-agent --agent my-agent --task-list experiments/splits/drug-discovery.txt

# All available tasks
biomlbench run-agent --agent my-agent --task-list experiments/splits/all.txt
```

## 🧪 Development

### Adding New Tasks

1. **Create task directory**:
```bash
mkdir -p biomlbench/tasks/my-new-task
```

2. **Implement required files**:
   - `config.yaml`: Task configuration
   - `prepare.py`: Data preparation logic
   - `grade.py`: Evaluation function
   - `description.md`: Task description
  
3. **Test the task**:
```bash
biomlbench prepare -t my-new-task
biomlbench run-agent --agent dummy --task-id my-new-task
biomlbench grade-sample runs/dummy/my-new-task/submission.csv my-new-task
```

### Building Custom Environments

To add specialized dependencies:

```dockerfile
FROM biomlbench-env

# Add your biomedical tools
RUN conda install -c conda-forge your-bio-package
RUN pip install specialized-ml-library

# Your agent code here
```

### Contributing

1. Fork the repository
2. Create a feature branch
3. Add your task or improvement
4. Test with the dummy agent
5. Submit a pull request

For questions or contributions, please open an issue on GitHub.
