from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Any, Optional

from appdirs import user_cache_dir

from biomlbench.grade_helpers import Grader
from biomlbench.utils import get_logger, get_module_dir, get_repo_dir, import_fn, load_yaml

logger = get_logger(__name__)


DEFAULT_DATA_DIR = (Path(user_cache_dir()) / "bioml-bench" / "data").resolve()


@dataclass(frozen=True)
class Task:
    """
    Represents a biomedical ML task in the BioML-bench framework.
    
    Extended from MLE-bench Competition class with biomedical-specific metadata.
    """
    id: str
    name: str
    description: str
    grader: Grader
    answers: Path
    gold_submission: Path
    sample_submission: Path
    task_type: str  # e.g., 'medical_imaging', 'protein_engineering'
    domain: str     # e.g., 'oncology', 'drug_discovery'
    difficulty: str # e.g., 'easy', 'medium', 'hard'
    prepare_fn: Callable[[Path, Path, Path], Path]
    raw_dir: Path
    private_dir: Path
    public_dir: Path
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
        assert isinstance(self.answers, Path), "Task answers must be a Path."
        assert isinstance(self.gold_submission, Path), "Gold submission must be a Path."
        assert isinstance(self.sample_submission, Path), "Sample submission must be a Path."
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

        try:
            return Task(
                id=data["id"],
                name=data["name"],
                description=data["description"],
                grader=grader,
                answers=data["answers"],
                sample_submission=data["sample_submission"],
                gold_submission=data["gold_submission"],
                task_type=data["task_type"],
                domain=data["domain"],
                difficulty=data["difficulty"],
                prepare_fn=data["prepare_fn"],
                raw_dir=data["raw_dir"],
                public_dir=data["public_dir"],
                private_dir=data["private_dir"],
                checksums=data["checksums"],
                leaderboard=data["leaderboard"],
                biomedical_metadata=data.get("biomedical_metadata"),
                human_baselines=data.get("human_baselines"),
                compute_requirements=data.get("compute_requirements"),
                data_source=data.get("data_source"),
            )
        except KeyError as e:
            raise ValueError(f"Missing key {e} in task config!")


class Registry:
    def __init__(self, data_dir: Path = DEFAULT_DATA_DIR):
        self._data_dir = data_dir.resolve()

    def get_task(self, task_id: str) -> Task:
        """Fetch the task from the registry."""

        config_path = self.get_tasks_dir() / task_id / "config.yaml"
        config = load_yaml(config_path)

        checksums_path = self.get_tasks_dir() / task_id / "checksums.yaml"
        leaderboard_path = self.get_tasks_dir() / task_id / "leaderboard.csv"

        description_path = get_repo_dir() / config["description"]
        description = description_path.read_text()

        preparer_fn = import_fn(config["preparer"])

        answers = self.get_data_dir() / config["dataset"]["answers"]
        gold_submission = answers
        if "gold_submission" in config["dataset"]:
            gold_submission = self.get_data_dir() / config["dataset"]["gold_submission"]
        sample_submission = self.get_data_dir() / config["dataset"]["sample_submission"]

        raw_dir = self.get_data_dir() / task_id / "raw"
        private_dir = self.get_data_dir() / task_id / "prepared" / "private"
        public_dir = self.get_data_dir() / task_id / "prepared" / "public"

        return Task.from_dict(
            {
                **config,
                "description": description,
                "answers": answers,
                "sample_submission": sample_submission,
                "gold_submission": gold_submission,
                "prepare_fn": preparer_fn,
                "raw_dir": raw_dir,
                "private_dir": private_dir,
                "public_dir": public_dir,
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
        task_ids = []
        for task_id in self.list_task_ids():
            task = self.get_task(task_id)
            if task.domain == domain:
                task_ids.append(task_id)
        return task_ids

    def get_tasks_by_type(self, task_type: str) -> list[str]:
        """List all task IDs for a specific task type."""
        task_ids = []
        for task_id in self.list_task_ids():
            task = self.get_task(task_id)
            if task.task_type == task_type:
                task_ids.append(task_id)
        return task_ids

    def get_lite_task_ids(self) -> list[str]:
        """List all task IDs for the lite version (low difficulty tasks)."""
        lite_tasks_file = self.get_splits_dir() / "lite.txt"
        if lite_tasks_file.exists():
            with open(lite_tasks_file, "r") as f:
                task_ids = f.read().splitlines()
            return task_ids
        else:
            # Fallback: return easy/medium difficulty tasks
            return self.get_tasks_by_difficulty("easy") + self.get_tasks_by_difficulty("medium")

    def get_tasks_by_difficulty(self, difficulty: str) -> list[str]:
        """List all task IDs for a specific difficulty level."""
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
        task_ids = [f.parent.stem for f in sorted(task_configs)]
        return task_ids




registry = Registry()
