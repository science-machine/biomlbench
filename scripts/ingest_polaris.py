#!/usr/bin/env python3
"""
Simple script to ingest a few Polaris tasks.
"""

import signal
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import TimeoutError as FuturesTimeoutError
from concurrent.futures import as_completed
from functools import partial
from pathlib import Path


class TaskTimeoutError(Exception):
    pass


def create_polaris_task(benchmark_id: str, output_dir: Path, prepare_data: bool = True):
    """Create a single Polaris task using PolarisDataSource."""

    # Generate task name
    task_name = benchmark_id.replace("/", "-").lower()
    task_dir = output_dir / task_name
    task_dir.mkdir(parents=True, exist_ok=True)

    print(f"Creating task: polarishub/{task_name}")

    # Use PolarisDataSource to get metadata without downloading
    from biomlbench.data_sources.polaris import PolarisDataSource

    data_source = PolarisDataSource()
    source_config = {"benchmark_id": benchmark_id}

    # Validate the config works
    data_source.validate_config(source_config)

    # Get actual metadata from Polaris
    import polaris as po

    benchmark = po.load_benchmark(benchmark_id)
    actual_name = benchmark.name
    actual_description = benchmark.description
    main_metric = benchmark.main_metric
    task_type_polaris = benchmark.task_type

    # Get data structure to understand targets
    train, _ = benchmark.get_train_test_split()
    train_df = train.as_dataframe()

    # First column is ALWAYS the feature/molecule column
    molecule_col = train_df.columns[0]

    # Use Polaris target_cols metadata to identify targets
    polaris_target_cols = benchmark.target_cols

    # Use first target if multiple targets exist
    if isinstance(polaris_target_cols, list) and len(polaris_target_cols) > 1:
        first_target = polaris_target_cols[0]
        print(f"  Multiple targets found: {polaris_target_cols}")
        print(f"  Using first target only: {first_target}")
    elif isinstance(polaris_target_cols, list) and len(polaris_target_cols) == 1:
        first_target = polaris_target_cols[0]
    else:
        # Single target (not a list)
        first_target = polaris_target_cols

    print(f"  All columns: {list(train_df.columns)}")
    print(f"  Molecule column: {molecule_col}")
    print(f"  Polaris target_cols: {polaris_target_cols}")
    print(f"  Using first target: {first_target}")
    print(f"  Main metric: {main_metric}")

    # Create config.yaml
    config = f"""id: polarishub/{task_name}
name: "{actual_name}"
task_type: drug_discovery
domain: molecular_properties
difficulty: medium
awards_medals: true
prizes: null
description: biomlbench/tasks/polarishub/{task_name}/description.md

data_source:
  type: polaris
  benchmark_id: {benchmark_id}

dataset:
  answers: polarishub/{task_name}/prepared/private/answers.csv
  sample_submission: polarishub/{task_name}/prepared/public/sample_submission.csv

grader:
  name: polaris-metric
  grade_fn: biomlbench.tasks.polarishub.{task_name}.grade:grade

preparer: biomlbench.tasks.polarishub.{task_name}.prepare:prepare

biomedical_metadata:
  modality: "molecular"
  data_type: "{task_type_polaris}"
  clinical_relevance: "drug_discovery"
  polaris_main_metric: "{main_metric}"

compute_requirements:
  recommended_gpu_memory_gb: 4
  estimated_runtime_minutes: 15
  max_dataset_size_gb: 2
"""

    # Create prepare.py that uses PolarisDataSource pattern
    prepare = f'''"""Auto-generated preparation for {benchmark_id}."""
from pathlib import Path
import pandas as pd

def prepare(raw: Path, public: Path, private: Path) -> None:
    """
    Prepare {benchmark_id} dataset from Polaris data downloaded by PolarisDataSource.
    
    Args:
        raw: Directory with polaris_train_data.csv and polaris_test_data.csv
        public: Directory for public data (train.csv, test_features.csv)  
        private: Directory for private data (answers.csv)
    """
    
    # Load the Polaris data files (downloaded by PolarisDataSource)
    train_df = pd.read_csv(raw / 'polaris_train_data.csv')
    test_df = pd.read_csv(raw / 'polaris_test_data.csv')
    
    print(f"Loaded Polaris data: {{len(train_df)}} train, {{len(test_df)}} test samples")
    print(f"Columns: {{list(train_df.columns)}}")
    
    # First column is ALWAYS the feature/molecule column
    molecule_col = train_df.columns[0]
    
    # Load benchmark to get target metadata
    import polaris as po
    benchmark = po.load_benchmark("{benchmark_id}")
    polaris_target_cols = list(benchmark.target_cols)
    
    # Use first target if multiple targets exist
    if len(polaris_target_cols) > 1:
        target_col = polaris_target_cols[0]
        print(f"Multiple targets found: {{polaris_target_cols}}")
        print(f"Using first target only: {{target_col}}")
    elif len(polaris_target_cols) == 1:
        target_col = polaris_target_cols[0]
    else:
        raise ValueError(f"No targets found in benchmark {benchmark_id}")
    
    print(f"Using molecule column: {{molecule_col}}")
    print(f"Using target column: {{target_col}}")
    
    # Create public training data (train.csv)
    train_public = train_df[[molecule_col, target_col]].copy()
    train_public.to_csv(public / 'train.csv', index=False)
    
    # Create test features (test_features.csv) - no targets
    test_features = test_df[[molecule_col]].copy()
    test_features['id'] = range(len(test_features))
    test_features = test_features[['id', molecule_col]]
    test_features.to_csv(public / 'test_features.csv', index=False)
    
    # Create sample submission (sample_submission.csv)
    sample_submission = pd.DataFrame({{
        'id': range(len(test_features)),
        target_col: [0.0] * len(test_features)
    }})
    sample_submission.to_csv(public / 'sample_submission.csv', index=False)
    
    # Create private answers file for evaluation
    answers = pd.DataFrame({{
        'id': range(len(test_features)),
        target_col: test_df[target_col].values
    }})
    answers.to_csv(private / 'answers.csv', index=False)
    
    print("✅ Dataset prepared successfully!")
    print(f"Molecule column: {{molecule_col}}")
    print(f"Target column: {{target_col}}")
    print(f"Files created:")
    print(f"  Public: train.csv, test_features.csv, sample_submission.csv")
    print(f"  Private: answers.csv")
'''

    # Create grade.py that uses main_metric
    grade = f'''"""Auto-generated grading for {benchmark_id}."""
import pandas as pd
import polaris as po

def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """Grade using Polaris evaluation with main_metric."""
    benchmark = po.load_benchmark("{benchmark_id}")
    
    # Extract predictions (second column)
    predictions = submission[submission.columns[1]].tolist()
    
    # Use Polaris evaluation
    try:
        results = benchmark.evaluate(predictions)
    except ValueError as e:
        if "Metric requires `y_prob` input" in str(e):
            results = benchmark.evaluate(y_prob=predictions)
        else:
            raise e
    
    # Return main metric score
    main_metric = "{main_metric}"
    
    if hasattr(results, 'results') and not results.results.empty:
        # Search for the main metric in the results
        if "Metric" in results.results.columns:
            # Find the row that matches the main metric
            metric_mask = results.results["Metric"] == main_metric
            matching_rows = results.results[metric_mask]
            
            if matching_rows.empty:
                raise ValueError(
                    f"Main metric '{main_metric}' not found in evaluation results. "
                )
            elif len(matching_rows) > 1:
                raise ValueError(
                    f"Multiple rows found for main metric '{main_metric}'. "
                )
            else:
                # Extract score from the single matching row
                score = float(matching_rows["Score"].iloc[0])
        else:
            raise ValueError(
                f"'Metric' column not found in evaluation results. "
            )
    else:
        raise ValueError(
            "Evaluation results are empty or missing. "
        )
    
    return score
'''

    # Create description.md with actual Polaris README
    polaris_readme = benchmark.readme or ""
    description = f"""# {actual_name}

{polaris_readme}

---

**Source:** [Polaris Hub - {benchmark_id}](https://polarishub.io)  
**Task Type:** {task_type_polaris}  
**Main Metric:** {main_metric}

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `{molecule_col}`
- Target column: `{first_target}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `{benchmark_id}`.
Main metric: **{main_metric}**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
"""

    # Write files
    (task_dir / "config.yaml").write_text(config)
    (task_dir / "prepare.py").write_text(prepare)
    (task_dir / "grade.py").write_text(grade)
    (task_dir / "description.md").write_text(description)
    (task_dir / "__init__.py").write_text("# Auto-generated Polaris task\n")

    print(f"✅ Created polarishub/{task_name}")

    # Optionally prepare the data immediately
    if prepare_data:
        try:
            from biomlbench.data import download_and_prepare_dataset
            from biomlbench.registry import registry

            print(f"🔄 Preparing data for polarishub/{task_name}...")
            task = registry.get_task(f"polarishub/{task_name}")
            download_and_prepare_dataset(task)
            print(f"✅ Prepared data for polarishub/{task_name}")
        except Exception as e:
            print(f"❌ Failed to prepare {task_name}: {e}")
            # Don't raise - let the task creation succeed even if preparation fails


def main():
    """Discover and ingest all Polaris tasks."""
    # Load number of workers from command line
    import argparse

    parser = argparse.ArgumentParser(description="Ingest Polaris tasks")
    parser.add_argument(
        "--workers", type=int, default=4, help="Number of workers to use for parallel processing"
    )
    parser.add_argument(
        "--prepare",
        action="store_true",
        default=True,
        help="Also prepare datasets after creating tasks (default: True)",
    )
    parser.add_argument(
        "--no-prepare",
        dest="prepare",
        action="store_false",
        help="Skip dataset preparation (only create task structure)",
    )
    args = parser.parse_args()

    output_dir = Path("biomlbench/tasks/polarishub")

    print("🔍 Discovering benchmarks from Polaris Hub...")

    # Discover all benchmarks using PolarisHubClient
    from polaris.hub.client import PolarisHubClient

    with PolarisHubClient() as client:
        benchmarks = []
        for chunk in range(10):
            result = client.list_benchmarks(offset=chunk * 100, limit=100)
            if len(result) == 0:
                break
            benchmarks.extend(result)

    # Only keep the benchmarks which are in scripts/polarishub_tasks.csv
    with open("scripts/polarishub_tasks.csv", "r") as f:
        polarishub_tasks = [line.strip() for line in f.readlines()]
        # Remove everything before the FIRST -
        # EVERYTHING AFTER THE FIRST - IS THE TASK NAME
        polarishub_tasks = [task.split("-", 1)[1] for task in polarishub_tasks]
    benchmarks = [benchmark for benchmark in benchmarks if any(task in benchmark for task in polarishub_tasks)]
    print(f"📊 Found {len(benchmarks)} benchmarks on Polaris Hub" )

    if args.prepare:
        print("🔄 Will also prepare datasets after creating tasks")
    else:
        print("⏭️  Will skip dataset preparation")

    # Create partial function to avoid pickling issues
    task_func = partial(create_polaris_task, output_dir=output_dir, prepare_data=args.prepare)

    # Use ProcessPoolExecutor for parallel processing
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        # Submit all tasks
        future_to_benchmark = {
            executor.submit(task_func, benchmark_id): benchmark_id for benchmark_id in benchmarks
        }

        # Process completed tasks with individual timeouts
        for future in as_completed(future_to_benchmark):
            benchmark_id = future_to_benchmark[future]
            try:
                # Set 10-minute timeout for each individual task
                future.result(timeout=600)  # 10 minutes per task
                status = "✅ Created and prepared" if args.prepare else "✅ Created"
                print(f"{status} polarishub/{benchmark_id.replace('/', '-').lower()}")
            except FuturesTimeoutError:
                print(f"⏰ Task timed out after 10 minutes: {benchmark_id}")
            except Exception as e:
                print(f"❌ Failed to create {benchmark_id}: {e}")


if __name__ == "__main__":
    main()
