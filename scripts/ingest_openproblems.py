#!/usr/bin/env python3
"""
Script to ingest OpenProblems batch_integration task.
Similar to ingest_polaris.py but for OpenProblems data.
"""

import sys
from pathlib import Path

# Add biomlbench to path
sys.path.insert(0, str(Path(__file__).parent.parent))

def create_openproblems_task(task_name: str, output_dir: Path, prepare_data: bool = True):
    """Create OpenProblems batch_integration task using OpenProblemsDataSource."""
    
    task_dir = output_dir / task_name
    task_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Creating task: openproblems/{task_name}")
    
    # Use OpenProblemsDataSource to get metadata
    from biomlbench.data_sources.openproblems import OpenProblemsDataSource
    
    data_source: OpenProblemsDataSource = OpenProblemsDataSource()
    source_config = {"task_name": "batch_integration"}
    
    # Download actual OpenProblems data to get metadata
    print("🔄 Downloading OpenProblems metadata from S3...")
    breakpoint()
    task_data = data_source.download_task_results("batch_integration")
    
    # Get task info
    task_info = task_data['task_info']
    actual_name = task_info['name']
    actual_description = task_info['summary']
    
    print(f"  Task name: {actual_name}")
    print(f"  Task summary: {actual_description}")
    print(f"  Available datasets: {len(task_data['datasets'])}")
    print(f"  Score entries: {len(task_data['scores'])}")
    
    # Create config.yaml
    config = f"""id: openproblems/{task_name}
name: "{actual_name}"
task_type: benchmark_results
domain: single_cell
difficulty: reference
awards_medals: false
prizes: null
description: biomlbench/tasks/openproblems/{task_name}/description.md

data_source:
  type: openproblems
  task_name: batch_integration

dataset:
  answers: openproblems/{task_name}/prepared/private/answers.csv
  sample_submission: openproblems/{task_name}/prepared/public/sample_submission.csv

grader:
  name: reference-only
  grade_fn: biomlbench.tasks.openproblems.{task_name}.grade:grade

preparer: biomlbench.tasks.openproblems.{task_name}.prepare:prepare

biomedical_metadata:
  modality: single_cell
  data_type: benchmark_results
  source: openproblems
  openproblems_task: batch_integration

compute_requirements:
  recommended_gpu_memory_gb: 0
  estimated_runtime_minutes: 1
  max_dataset_size_gb: 0.001
"""

    # Create prepare.py that demonstrates the data flow
    prepare = f'''"""OpenProblems batch_integration data preparation."""
from pathlib import Path
import pandas as pd
import logging

logger = logging.getLogger(__name__)

def prepare(raw: Path, public: Path, private: Path) -> None:
    """
    Prepare OpenProblems batch_integration dataset.
    
    This is a results-only task that displays pre-computed benchmark results.
    The "data preparation" creates placeholder files to satisfy biomlbench structure.
    
    Args:
        raw: Directory for raw data (unused for results-only tasks)
        public: Directory for public data (sample submission, description)
        private: Directory for private data (answers)
    """
    
    logger.info("This is a reference-only task from OpenProblems")
    logger.info("Creating placeholder files to satisfy biomlbench structure...")
    
    # Download the actual results to show in logs
    from biomlbench.data_sources.openproblems import OpenProblemsDataSource
    
    print("🔄 Downloading OpenProblems results from S3...")
    data_source = OpenProblemsDataSource()
    task_data = data_source.download_task_results("batch_integration")
    
    # Get the leaderboard
    leaderboard = data_source.create_simple_leaderboard(task_data)
    detailed_results = data_source.get_detailed_results({{"task_name": "batch_integration"}})
    
    print(f"📊 Downloaded {{len(detailed_results)}} method-metric-dataset combinations")
    print(f"🏆 Leaderboard has {{len(leaderboard)}} methods")
    print("\\nTop methods:")
    for _, row in leaderboard.head(3).iterrows():
        print(f"  {{row['teamName']}}: {{row['score']:.4f}}")
    
    # Create sample submission (required by biomlbench structure)
    sample_submission = pd.DataFrame({{
        'id': [0, 1, 2],
        'prediction': [0.5, 0.6, 0.7]
    }})
    sample_submission.to_csv(public / 'sample_submission.csv', index=False)
    
    # Create answers file (required by biomlbench structure)  
    answers = pd.DataFrame({{
        'id': [0, 1, 2],
        'prediction': [0.5, 0.6, 0.7]
    }})
    answers.to_csv(private / 'answers.csv', index=False)
    
    # Create task description for public access
    description_content = f"""# OpenProblems Batch Integration Results

This task displays pre-computed benchmark results from OpenProblems.

## Leaderboard

{{leaderboard.to_markdown(index=False)}}

## About

- **Data source**: OpenProblems S3 bucket
- **Task type**: Results-only benchmark
- **Total method-metric combinations**: {{len(detailed_results)}}
- **Unique methods**: {{detailed_results['method'].nunique()}}
- **Unique metrics**: {{detailed_results['metric'].nunique()}}
- **Unique datasets**: {{detailed_results['dataset'].nunique()}}

Visit https://openproblems.bio for actual evaluation and submission.
"""
    
    (public / 'description.md').write_text(description_content)
    
    logger.info("✅ Placeholder files created successfully")
    logger.info("This task displays OpenProblems benchmark results")
    logger.info("Visit https://openproblems.bio for actual evaluation")
'''

    # Create grade.py for results-only grading
    grade = f'''"""OpenProblems batch_integration grading (reference-only)."""
import pandas as pd
import logging

logger = logging.getLogger(__name__)

def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
    """
    Reference-only grading for OpenProblems batch_integration.
    
    This task displays pre-computed results and doesn't perform actual grading.
    """
    
    logger.info("This is a reference-only task from OpenProblems")
    logger.info("No grading available - results are pre-computed")
    logger.info("Visit https://openproblems.bio for actual evaluation")
    
    # Show the current leaderboard when "grading" is attempted
    from biomlbench.data_sources.openproblems import OpenProblemsDataSource
    
    data_source = OpenProblemsDataSource()
    leaderboard = data_source.get_leaderboard({{"task_name": "batch_integration"}})
    
    print("\\n🏆 Current OpenProblems Leaderboard:")
    print(leaderboard.to_string(index=False))
    
    # Return 0 for reference-only tasks
    return 0.0
'''

    # Create description.md with actual OpenProblems info
    task_info = task_data['task_info']
    description_md = f"""# {actual_name}

## Summary
{task_info['summary']}

## Description
{task_info.get('description', 'Remove unwanted batch effects from scRNA data while retaining biologically meaningful variation.')}

## Motivation
{task_info.get('motivation', 'Batch effects are technical artifacts that can confound biological signal in single-cell data analysis.')}

---

**Source:** [OpenProblems - {task_name}](https://openproblems.bio)  
**Task Type:** benchmark_results  
**Data Type:** Results-only (pre-computed)

## Available Data

This task provides access to pre-computed benchmark results from OpenProblems:

- **Methods**: {len(set(entry['method_id'] for entry in task_data['scores']))} integration methods
- **Datasets**: {len(task_data['datasets'])} single-cell datasets
- **Metrics**: {len(set(metric for entry in task_data['scores'] for metric in entry['metric_ids']))} evaluation metrics
- **Total Results**: {len(task_data['scores'])} method-dataset combinations

## Evaluation Metrics

OpenProblems uses multiple metrics to evaluate batch integration:
- **KBET**: k-Nearest Neighbor Batch Effect Test
- **ASW (batch)**: Average Silhouette Width for batch mixing
- **ASW (label)**: Average Silhouette Width for label separation  
- **Graph connectivity**: Preservation of local neighborhoods
- **ARI**: Adjusted Rand Index
- **NMI**: Normalized Mutual Information
- **Isolated label ASW**: Silhouette width for isolated labels
- **Isolated label F1**: F1 score for isolated label detection

## Usage

This is a **results-only** task that displays pre-computed benchmarks. For actual submission and evaluation, visit [OpenProblems](https://openproblems.bio).

## Source

Auto-generated from [OpenProblems](https://openproblems.bio/).
"""

    # Write all files
    (task_dir / "config.yaml").write_text(config)
    (task_dir / "prepare.py").write_text(prepare)
    (task_dir / "grade.py").write_text(grade)
    (task_dir / "description.md").write_text(description_md)
    (task_dir / "__init__.py").write_text("# Auto-generated OpenProblems task\\n")
    
    # Save the raw results data for reference
    import yaml
    (task_dir / "results_data.yaml").write_text(yaml.dump(task_data, default_flow_style=False))
    
    # Create the leaderboard file
    breakpoint()
    leaderboard = data_source.create_simple_leaderboard(task_data)
    leaderboard.to_csv(task_dir / "leaderboard.csv", index=False)
    
    print(f"✅ Created openproblems/{task_name}")
    print(f"📊 Leaderboard: {len(leaderboard)} methods")
    print(f"📁 Results data: {len(task_data['scores'])} entries")
    
    # Optionally prepare the data immediately (run the prepare pipeline)
    if prepare_data:
        try:
            from biomlbench.data import download_and_prepare_dataset
            from biomlbench.registry import registry
            
            print(f"🔄 Running data preparation pipeline for openproblems/{task_name}...")
            task = registry.get_task(f"openproblems/{task_name}")
            download_and_prepare_dataset(task)
            print(f"✅ Data preparation completed for openproblems/{task_name}")
            
            # Show what was created  
            print(f"📁 Data directory: ~/.cache/bioml-bench/data/openproblems/{task_name}")
            
            # List created files in the cache directory
            import os
            cache_dir = Path.home() / ".cache" / "bioml-bench" / "data" / "openproblems" / task_name
            if cache_dir.exists():
                print("📄 Created files:")
                for file_path in sorted(cache_dir.rglob("*")):
                    if file_path.is_file():
                        rel_path = file_path.relative_to(cache_dir)
                        print(f"  {rel_path}")
            
        except Exception as e:
            print(f"❌ Failed to prepare {task_name}: {e}")
            import traceback
            traceback.print_exc()
            # Don't raise - let the task creation succeed even if preparation fails


def main():
    """Create and prepare the OpenProblems batch_integration task."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Ingest OpenProblems batch_integration task")
    parser.add_argument(
        "--prepare",
        action="store_true", 
        default=True,
        help="Also run data preparation pipeline after creating task (default: True)"
    )
    parser.add_argument(
        "--no-prepare", 
        dest="prepare",
        action="store_false",
        help="Skip data preparation (only create task structure)"
    )
    args = parser.parse_args()
    
    output_dir = Path("biomlbench/tasks/openproblems")
    task_name = "batch_integration"
    
    print("🔬 Creating OpenProblems batch_integration task...")
    
    if args.prepare:
        print("🔄 Will also run data preparation pipeline")
    else:
        print("⏭️  Will skip data preparation pipeline")
    
    # Create the task
    create_openproblems_task(task_name, output_dir, args.prepare)
    
    print("\\n✅ OpenProblems task creation complete!")
    print(f"📁 Task location: biomlbench/tasks/openproblems/{task_name}")
    print("🔬 This demonstrates the full OpenProblems data pipeline")


if __name__ == "__main__":
    main() 