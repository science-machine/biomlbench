from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Optional

from appdirs import user_cache_dir

from biomlbench.grade_helpers import Grader
from biomlbench.utils import (
    get_logger,
    get_module_dir,
    get_repo_dir,
    import_fn,
    is_empty,
    load_yaml,
)

logger = get_logger(__name__)


DEFAULT_DATA_DIR = (Path(user_cache_dir()) / "bioml-bench" / "data").resolve()


@dataclass(frozen=True)
class Dataset:
    """
    Represents a dataset in the BioML-bench framework.

    Extended from MLE-bench CompetitionDataset class with biomedical-specific metadata.
    """

    task_id: str  # ID of the task to which the dataset belongs
    dataset_id: str  # ID of the dataset
    path: Path  # path to the dataset directory containing public and private subdirectories
    answers: str  # name of the answers file
    sample_submission: str  # name of the sample submission file

    @property
    def id(self) -> str:
        """Combines the task ID and the dataset ID to form a unique identifier for the dataset."""
        return f"{self.task_id}-{self.dataset_id}"

    @property
    def answers_path(self) -> Path:
        """Path to the answers file."""
        return self.path / "private" / self.answers

    @property
    def public_path(self) -> Path:
        """Path to the public directory."""
        return self.path / "public"

    @property
    def private_path(self) -> Path:
        """Path to the private directory."""
        return self.path / "private"

    @property
    def sample_submission_path(self) -> Path:
        """Path to the sample submission file."""
        return self.path / "public" / self.sample_submission

    def is_prepared(self) -> bool:
        """Checks if the dataset is prepared by checking that the public and private directories are not empty."""
        return not (is_empty(self.public_path) or is_empty(self.private_path))


@dataclass(frozen=True)
class Task:
    """
    Represents a biomedical ML task in the BioML-bench framework, which can contain multiple datasets.

    Extended from MLE-bench Competition class with biomedical-specific metadata.
    """

    id: str
    name: str
    description: str
    datasets: list[Dataset]
    grader: Grader
    task_type: str  # e.g., 'medical_imaging', 'protein_engineering'
    domain: str  # e.g., 'oncology', 'drug_discovery'
    difficulty: str  # e.g., 'easy', 'medium', 'hard'
    prepare_fn: Callable[[Path, Path], None]
    raw_dir: Path
    prepared_dir: Path
    checksums: Path
    leaderboard: Path

    # Biomedical-specific fields
    biomedical_metadata: Optional[Dict[str, Any]] = None
    human_baselines: Optional[Dict[str, Any]] = None
    compute_requirements: Optional[Dict[str, Any]] = None
    data_source: Optional[Dict[str, Any]] = None

    def __post_init__(self):
        assert isinstance(self.id, str), "Task id must be a string."
        assert isinstance(self.name, str), "Task name must be a string."
        assert isinstance(self.description, str), "Task description must be a string."
        assert isinstance(self.grader, Grader), "Task grader must be of type Grader."
        assert isinstance(self.task_type, str), "Task type must be a string."
        assert isinstance(self.domain, str), "Domain must be a string."
        assert isinstance(self.difficulty, str), "Difficulty must be a string."
        assert isinstance(self.checksums, Path), "Checksums must be a Path."
        assert isinstance(self.leaderboard, Path), "Leaderboard must be a Path."
        assert len(self.id) > 0, "Task id cannot be empty."
        assert len(self.name) > 0, "Task name cannot be empty."
        assert len(self.description) > 0, "Task description cannot be empty."
        assert len(self.task_type) > 0, "Task type cannot be empty."
        assert len(self.domain) > 0, "Domain cannot be empty."
        assert len(self.difficulty) > 0, "Difficulty cannot be empty."

    @staticmethod
    def from_dict(data: dict) -> "Task":
        grader = Grader.from_dict(data["grader"])

        # Create dataset objects
        datasets = []

        # dataset_data here is expected to contain answers and sample_submission file names
        for dataset_id, dataset_data in data["datasets"].items():
            dataset = Dataset(
                task_id=data["id"],
                dataset_id=dataset_id,
                path=data["prepared_dir"] / dataset_id,
                answers=dataset_data["answers"],
                sample_submission=dataset_data["sample_submission"],
            )
            datasets.append(dataset)

        try:
            if not "task_type" in data:
                data["task_type"] = "simple"
            if not "domain" in data:
                data["domain"] = "general"
            if not "difficulty" in data:
                data["difficulty"] = "easy"

            return Task(
                id=data["id"],
                name=data["name"],
                description=data["description"],
                datasets=datasets,
                grader=grader,
                task_type=data["task_type"],
                domain=data["domain"],
                difficulty=data["difficulty"],
                prepare_fn=data["prepare_fn"],
                raw_dir=data["raw_dir"],
                prepared_dir=data["prepared_dir"],
                checksums=data["checksums"],
                leaderboard=data["leaderboard"],
                biomedical_metadata=data.get("biomedical_metadata"),
                human_baselines=data.get("human_baselines"),
                compute_requirements=data.get("compute_requirements"),
                data_source=data.get("data_source"),
            )
        except KeyError as e:
            raise ValueError(f"Missing key {e} in task config!")

    def get_dataset(self, dataset_id: str) -> Dataset:
        """Retrieves the dataset from the task using folder/dataset format."""
        for dataset in self.datasets:
            if dataset.id == dataset_id:
                return dataset
        raise ValueError(f"Dataset {dataset_id} not found in task {self.id}")

    def is_prepared(self) -> bool:
        """Checks if the task is prepared by checking that each of its datasets is prepared."""
        for dataset in self.datasets:
            try:
                if not dataset.is_prepared():
                    return False
            except AssertionError:
                return False

        return True


class Registry:
    def __init__(self, data_dir: Path = DEFAULT_DATA_DIR):
        self._data_dir = data_dir.resolve()

    def get_task(self, task_id: str) -> Task:
        """Fetch the task from the registry using folder/task format."""

        # New format: folder/task
        try:
            folder, task_name = task_id.split("/", 1)
        except ValueError:
            raise ValueError(f"Invalid task ID: {task_id}. Expected format: folder/task")

        config_path = self.get_tasks_dir() / folder / task_name / "config.yaml"
        checksums_path = self.get_tasks_dir() / folder / task_name / "checksums.yaml"
        leaderboard_path = self.get_tasks_dir() / folder / task_name / "leaderboard.csv"
        description_path = self.get_tasks_dir() / folder / task_name / "description.md"

        config = load_yaml(config_path)
        description = description_path.read_text()

        preparer_fn = import_fn(config["preparer"])

        raw_dir = self.get_data_dir() / task_id / "raw"
        prepared_dir = self.get_data_dir() / task_id / "prepared"

        # {dataset_id: {answers: str, sample_submission: str}}
        datasets = config["datasets"]

        return Task.from_dict(
            {
                **config,
                "datasets": datasets,
                "description": description,
                "prepare_fn": preparer_fn,
                "raw_dir": raw_dir,
                "prepared_dir": prepared_dir,
                "checksums": checksums_path,
                "leaderboard": leaderboard_path,
            }
        )

    def get_tasks_dir(self) -> Path:
        """Retrieves the tasks directory within the registry."""
        return get_module_dir() / "tasks"

    def get_splits_dir(self) -> Path:
        """Retrieves the splits directory within the repository."""
        return get_repo_dir() / "experiments" / "splits"

    def get_tasks_by_domain(self, domain: str) -> list[str]:
        """List all task IDs for a specific biomedical domain."""
        raise NotImplementedError("Domain filtering not implemented yet.")

    def get_tasks_by_type(self, task_type: str) -> list[str]:
        """List all task IDs for a specific task type."""
        raise NotImplementedError("Task type filtering not implemented yet.")

    def get_lite_task_ids(self) -> list[str]:
        """List all task IDs for the lite version (low difficulty tasks)."""
        raise NotImplementedError("Lite version not implemented yet.")
        lite_tasks_file = self.get_splits_dir() / "lite.txt"
        if lite_tasks_file.exists():
            with open(lite_tasks_file, "r") as f:
                task_ids = f.read().splitlines()
            return task_ids
        else:
            # Fallback: return easy/medium difficulty tasks
            raise NotImplementedError("Lite version not implemented yet.")

    def get_tasks_by_difficulty(self, difficulty: str) -> list[str]:
        """List all task IDs for a specific difficulty level."""
        raise NotImplementedError("Difficulty filtering not implemented yet.")
        task_ids = []
        for task_id in self.list_task_ids():
            task = self.get_task(task_id)
            if task.difficulty == difficulty:
                task_ids.append(task_id)
        return task_ids

    def get_data_dir(self) -> Path:
        """Retrieves the data directory within the registry."""
        return self._data_dir

    def set_data_dir(self, new_data_dir: Path) -> "Registry":
        """Sets the data directory within the registry."""
        return Registry(new_data_dir)

    def list_task_ids(self) -> list[str]:
        """List all task IDs available in the registry, sorted alphabetically."""
        task_configs = self.get_tasks_dir().rglob("config.yaml")
        task_ids = []

        for config_path in sorted(task_configs):
            # Check relative path from tasks directory
            rel_path = config_path.relative_to(self.get_tasks_dir())
            parts = rel_path.parts

            if len(parts) == 3:  # folder/task/config.yaml
                folder_name = parts[0]
                task_name = parts[1]
                if folder_name not in ["__pycache__"]:
                    task_ids.append(f"{folder_name}/{task_name}")
            elif len(parts) == 2:  # task/config.yaml (legacy)
                task_name = parts[0]
                if task_name not in ["__pycache__"]:
                    task_ids.append(task_name)

        return task_ids


registry = Registry()
