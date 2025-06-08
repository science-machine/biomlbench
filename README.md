# BioML-bench

**THIS WAS GENERATED IN COLLABORATION WITH CLAUDE. EXPECT BUGS.**

**A benchmark suite for evaluating machine learning agents on biomedical tasks**

BioML-bench is built on top of [MLE-bench](https://github.com/openai/mle-bench) and provides a comprehensive framework for benchmarking AI agents on biomedical machine learning (BioML) tasks including protein engineering, drug discovery, medical imaging, and clinical biomarkers.

~~Development TODOs (unordered)~~

- [ ] Add `experiments/splits` which will contain files like `medical-imaging.txt` that contain a list of task ids for only a given domain, or type of benchmarking run (e.g., lite, all).
- [ ] Decide if we want to have `experiments/familiarity` like the original MLE-bench. This would involve measuring the relationship between model performance and the "familiarity" of the model for the task description. This tries to assess if the model is "cheating" by already understanding the task.
- [ ] Need to decide if we can add human baselines in for Polaris. We would probably have to manually scrape the leaderboard (possibly with playwright) and then we would need to add a `supports_human_baselines` method to the `PolarisDataSource` class which returns `True`. And we would need to add a `get_human_baselines` method to the `PolarisDataSource` class which returns a pandas dataframe with the human baselines.
- [ ] Need an automated pipeline for scraping all the polaris and kaggle tasks and creating tasks/ directories for them.
- [ ] 


## How to wrap a new benchmark database

Example: Polaris.

To wrap polaris, we needed to:

1. Understand the data source. In this case, Polaris has a high-level API that can be used to automatically download and split data.
2. Create a new data source module in `biomlbench/data_sources/`. In the case of Polaris, we created `biomlbench/data_sources/polaris.py`.
3. Implement the `XYZDataSource` class. This class must implement the `download` and `get_leaderboard` methods.
4. 



## 🧬 Features

- **Diverse Biomedical Tasks**: Protein engineering, drug discovery, medical imaging, clinical biomarkers
- **Agent-Agnostic Evaluation**: Any agent that can produce CSV outputs can be evaluated
- **Standardized Metrics**: Domain-specific evaluation metrics (RMSD, TM-score, AUC-ROC, etc.)
- **Human Baselines**: Built-in human performance benchmarks for comparison
- **Secure Evaluation**: Containerized execution with no data leakage
- **Extensible Framework**: Easy to add new biomedical tasks

## 🚀 Quick Start

### Installation

BioML-bench uses [uv](https://github.com/astral-sh/uv) for dependency management:

```bash
# Install uv if you haven't already
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone the repository
git clone https://github.com/science-machine-learning/biomlbench.git
cd biomlbench

# Install dependencies
uv sync
```

### Basic Usage

1. **Prepare a task dataset:**
```bash
biomlbench prepare -t caco2-wang
```

2. **Run baselines to generate submissions:**
```bash
biomlbench run-baseline --baseline all caco2-wang
```

3. **Run evaluation on multiple results:**
```bash
biomlbench grade --submission submissions.jsonl --output-dir results/
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
- `caco2-wang`: Molecular property prediction, ADMET prediction
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

### For Agent Developers

To evaluate your agent on BioML-bench:

1. **Create a Docker container** that can read task data and produce submissions:
```dockerfile
FROM biomlbench-env
COPY my_agent.py /home/agent/
CMD ["python", "/home/agent/my_agent.py"]
```

2. **Agent Interface**: Your agent should:
   - Read task description from `/home/data/description.md`
   - Load training data from `/home/data/train.csv`
   - Load test features from `/home/data/test_features.csv`
   - Output predictions to `/home/submission/submission.csv`

3. **Run evaluation**:
```bash
python run_agent.py --agent-id my-agent --task-set experiments/splits/medical-imaging.txt
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
  
3. Optionally, implement baselines:
```bash
touch biomlbench/tasks/my-new-task/baselines.py
# Add your baselines to this file
```

4. **Test the task**:
```bash
biomlbench prepare -t my-new-task
biomlbench run-baseline --baseline all my-new-task
biomlbench grade-sample sample_submission.csv my-new-task
```
