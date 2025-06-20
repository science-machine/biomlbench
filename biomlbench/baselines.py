"""
Generic baseline framework for BioML-bench tasks.

This module provides the abstract baseline infrastructure. Task-specific
baselines should be implemented in each task directory (e.g., 
biomlbench/tasks/caco2-wang/baselines.py) and registered through the
baseline factory.
"""

import importlib.util
import json
import sys
import uuid
from abc import ABC, abstractmethod
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Type

import numpy as np
import pandas as pd

from biomlbench.data import is_dataset_prepared
from biomlbench.registry import Task
from biomlbench.utils import get_logger, get_module_dir, import_fn

logger = get_logger(__name__)


class BaselineAgent(ABC):
    """Abstract base class for baseline agents."""

    def __init__(self, name: str, seed: int = 42):
        self.name = name
        self.seed = seed
        np.random.seed(seed)

    @abstractmethod
    def generate_submission(self, task: Task) -> pd.DataFrame:
        """Generate a submission for the given task."""
        pass

    def get_metadata(self) -> Dict[str, Any]:
        """Get metadata about this baseline."""
        return {
            "baseline_name": self.name,
            "seed": self.seed,
            "timestamp": datetime.now().isoformat(),
        }


class SimpleBaselineAgent(BaselineAgent):
    """
    Generic simple baseline that predicts training data mean.
    Works for any regression task with numerical targets.
    """

    def __init__(self, seed: int = 42):
        super().__init__("simple", seed)

    def generate_submission(self, task: Task) -> pd.DataFrame:
        """Predict the mean of training data for all samples."""
        # Load sample submission to get format
        sample_submission = pd.read_csv(task.sample_submission)

        # Load training data
        train_df = pd.read_csv(task.public_dir / "train.csv")

        # Get target column (assume last column that's not 'id')
        # TODO: Revisit this logic as we add more data sources.
        target_cols = [col for col in train_df.columns if col.lower() != "id"]
        if len(target_cols) < 2:  # Need at least features + target
            raise ValueError("Could not identify target column in training data")

        # Assume last column is target (common convention)
        target_col = target_cols[-1]

        if not pd.api.types.is_numeric_dtype(train_df[target_col]):
            raise ValueError(f"Target column '{target_col}' is not numeric")

        mean_value = train_df[target_col].mean()

        # Create submission
        submission = sample_submission.copy()
        prediction_cols = [col for col in submission.columns if col.lower() != "id"]
        if not prediction_cols:
            raise ValueError("Could not identify prediction column in sample submission")

        prediction_col = prediction_cols[0]  # Use first prediction column
        submission[prediction_col] = mean_value

        logger.info(
            f"Generated {len(submission)} simple baseline predictions (mean={mean_value:.3f})"
        )
        return submission


class RandomBaselineAgent(BaselineAgent):
    """
    Generic random baseline that generates predictions based on training statistics.
    Works for any regression task with numerical targets.
    """

    def __init__(self, seed: int = 42):
        super().__init__("random", seed)

    def generate_submission(self, task: Task) -> pd.DataFrame:
        """Generate random predictions based on training data statistics."""
        # Load sample submission to get format
        sample_submission = pd.read_csv(task.sample_submission)

        # Load training data to get statistics
        train_df = pd.read_csv(task.public_dir / "train.csv")

        # Get target column (assume last column that's not 'id')
        target_cols = [col for col in train_df.columns if col.lower() != "id"]
        if len(target_cols) < 2:
            raise ValueError("Could not identify target column in training data")

        target_col = target_cols[-1]

        if not pd.api.types.is_numeric_dtype(train_df[target_col]):
            raise ValueError(f"Target column '{target_col}' is not numeric")

        target_stats = train_df[target_col].describe()

        # Generate random predictions within reasonable bounds
        n_predictions = len(sample_submission)
        mean = target_stats["mean"]
        std = target_stats["std"]

        # Generate predictions with some noise
        predictions = np.random.normal(mean, std * 0.5, n_predictions)

        # Create submission DataFrame
        submission = sample_submission.copy()
        prediction_cols = [col for col in submission.columns if col.lower() != "id"]
        if not prediction_cols:
            raise ValueError("Could not identify prediction column in sample submission")

        prediction_col = prediction_cols[0]
        submission[prediction_col] = predictions

        logger.info(
            f"Generated {n_predictions} random baseline predictions (mean={mean:.3f}, std={std:.3f})"
        )
        return submission


class BaselineFactory:
    """Factory for discovering and creating task-specific baselines."""

    _registry: Dict[str, Dict[str, Type[BaselineAgent]]] = {}

    @classmethod
    def register_baseline(
        cls, task_id: str, baseline_name: str, baseline_class: Type[BaselineAgent]
    ):
        """Register a baseline for a specific task."""
        if task_id not in cls._registry:
            cls._registry[task_id] = {}
        cls._registry[task_id][baseline_name] = baseline_class
        logger.debug(f"Registered {baseline_name} baseline for task {task_id}")

    @classmethod
    def get_available_baselines(cls, task_id: str) -> List[str]:
        """Get list of available baselines for a task."""
        # Always include generic baselines
        baselines = ["simple", "random"]

        # Add task-specific baselines if available
        if task_id in cls._registry:
            baselines.extend(cls._registry[task_id].keys())

        return list(set(baselines))  # Remove duplicates

    @classmethod
    def create_baseline(
        cls, task_id: str, baseline_type: str, seed: int = 42, **kwargs
    ) -> BaselineAgent:
        """Create a baseline agent for the given task and type."""
        # Try task-specific baselines first
        if task_id in cls._registry and baseline_type in cls._registry[task_id]:
            baseline_class = cls._registry[task_id][baseline_type]
            return baseline_class(seed=seed, **kwargs)

        # Fall back to generic baselines
        if baseline_type == "simple":
            return SimpleBaselineAgent(seed=seed)
        elif baseline_type == "random":
            return RandomBaselineAgent(seed=seed)
        else:
            available = cls.get_available_baselines(task_id)
            raise ValueError(
                f"Unknown baseline type '{baseline_type}' for task '{task_id}'. Available: {available}"
            )

    @classmethod
    def auto_discover_baselines(cls, task_id: str):
        """Auto-discover baselines for a task by trying to import task-specific baseline module."""
        try:
            # Construct the path to the task's baselines module
            task_dir = get_module_dir() / "tasks" / task_id
            baselines_file = task_dir / "baselines.py"

            if baselines_file.exists():
                # Create a module spec and load the module
                spec = importlib.util.spec_from_file_location(
                    f"task_{task_id.replace('/', '_').replace('-', '_')}_baselines", baselines_file
                )
                baseline_module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(baseline_module)

                # If the module has a register_baselines function, call it
                if hasattr(baseline_module, "register_baselines"):
                    baseline_module.register_baselines()
                    logger.info(f"Auto-discovered baselines for task '{task_id}'")
                else:
                    logger.debug(f"No register_baselines function found for task '{task_id}'")
            else:
                logger.debug(f"No baselines.py file found for task '{task_id}'")

        except Exception as e:
            logger.debug(f"Failed to auto-discover baselines for '{task_id}': {e}")


def run_baseline(
    task: Task,
    baseline_type: str = "simple",
    output_dir: Path = Path("baseline_submissions"),
    seed: int = 42,
    **kwargs,
) -> List[Dict[str, Any]]:
    """
    Run baseline(s) on a task and save submissions.

    Args:
        task: The task to run baselines on
        baseline_type: Type of baseline ('simple', 'random', or task-specific types, 'all')
        output_dir: Directory to save submissions
        seed: Random seed for reproducibility
        **kwargs: Additional arguments for specific baselines

    Returns:
        List of submission metadata dictionaries
    """
    # Verify task is prepared
    if not is_dataset_prepared(task):
        raise ValueError(
            f"Task '{task.id}' is not prepared. Run 'biomlbench prepare -t {task.id}' first."
        )

    # Auto-discover task-specific baselines
    BaselineFactory.auto_discover_baselines(task.id)

    output_dir.mkdir(parents=True, exist_ok=True)
    submission_metadata = []

    # Determine which baselines to run
    if baseline_type == "all":
        baseline_types = BaselineFactory.get_available_baselines(task.id)
    else:
        baseline_types = [baseline_type]

    for btype in baseline_types:
        logger.info(f"Running {btype} baseline for task '{task.id}'...")

        try:
            # Create baseline agent
            baseline = BaselineFactory.create_baseline(task.id, btype, seed=seed, **kwargs)

            # Generate submission
            submission = baseline.generate_submission(task)

            # Save submission
            run_id = str(uuid.uuid4())[:8]
            task_id_safe = task.id.replace("/", "-")
            filename = f"{task_id_safe}_{btype}_baseline_{run_id}.csv"
            submission_path = output_dir / filename
            submission.to_csv(submission_path, index=False)

            # Save metadata
            metadata = {
                "task_id": task.id,
                "baseline_type": btype,
                "submission_path": str(submission_path),
                "run_id": run_id,
                "seed": seed,
                **baseline.get_metadata(),
            }

            metadata_path = output_dir / f"{task_id_safe}_{btype}_baseline_{run_id}_metadata.json"
            with open(metadata_path, "w") as f:
                json.dump(metadata, f, indent=2)

            submission_metadata.append(metadata)

            logger.info(f"✅ {btype} baseline completed. Submission saved to: {submission_path}")

        except Exception as e:
            logger.error(f"❌ Failed to run {btype} baseline: {e}")
            continue

    # Create a summary JSONL file for easy grading
    if submission_metadata:
        task_id_safe = task.id.replace("/", "-")
        jsonl_path = output_dir / f"{task_id_safe}_baselines.jsonl"
        with open(jsonl_path, "w") as f:
            for metadata in submission_metadata:
                f.write(
                    json.dumps(
                        {
                            "task_id": metadata["task_id"],
                            "submission_path": metadata["submission_path"],
                        }
                    )
                    + "\n"
                )

        logger.info(f"📊 Summary file created: {jsonl_path}")
        logger.info(
            f"🎯 To grade all baselines, run: biomlbench grade --submission {jsonl_path} --output-dir {output_dir}"
        )

    return submission_metadata
