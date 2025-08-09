"""High-level grading functionality"""
import json
from datetime import datetime
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from biomlbench.data import get_leaderboard, is_dataset_prepared
from biomlbench.grade_helpers import TaskReport
from biomlbench.registry import Registry, Task
from biomlbench.registry import registry as DEFAULT_REGISTRY
from biomlbench.utils import get_logger, get_timestamp, load_answers, purple, read_csv, read_jsonl

logger = get_logger(__name__)


def grade_jsonl(
    path_to_submissions: Path,
    output_dir: Path,
    registry: Registry = DEFAULT_REGISTRY,
):
    """
    Grades multiple submissions stored in a JSONL file.
    Saves the aggregated report as a JSON file.
    """

    submissions = read_jsonl(str(path_to_submissions), skip_commented_out_lines=True)
    task_reports = []

    for submission in tqdm(submissions, desc="Grading submissions", unit="submission"):
        # Resolve submission path relative to the JSONL file directory
        submission_path = path_to_submissions.parent / submission["submission_path"]
        task_id = submission["task_id"]
        task = registry.get_task(task_id)
        single_report = grade_submission(submission_path, task)
        task_reports.append(single_report)

    aggregated_report = aggregate_reports(task_reports)
    timestamp = get_timestamp()
    save_path = output_dir / f"{timestamp}_grading_report.json"
    logger.info(
        json.dumps({k: v for k, v in aggregated_report.items() if k != "task_reports"}, indent=4)
    )

    output_dir.mkdir(exist_ok=True)
    with open(save_path, "w") as f:
        json.dump(aggregated_report, f, indent=2)
    logger.info(purple(f"Saved summary report to {save_path}"))


def grade_submission(path_to_submission: Path, task: Task) -> TaskReport:
    """Grades a submission file (CSV or h5ad) for the given task."""

    if not is_dataset_prepared(task, grading_only=True):
        raise ValueError(
            f"Dataset for task `{task.id}` is not prepared! "
            f"Please run `biomlbench prepare -t {task.id}` to prepare the dataset."
        )

    score = None
    submission_exists = path_to_submission.is_file()
    file_extension = path_to_submission.suffix.lower()

    if submission_exists:
        if file_extension == ".csv":
            # Traditional CSV-based task
            submission_df = read_csv(path_to_submission)
            answers = load_answers(task.answers)
            score = task.grader(submission_df, answers)
        elif file_extension == ".h5ad":
            # AnnData-based task (e.g., OpenProblems)
            # Call grade_fn directly with Path arguments, bypassing Grader.__call__
            answers_path = Path(task.answers)
            try:
                score = task.grader.grade_fn(path_to_submission, answers_path)
                if score is not None:
                    score = round(score, 5)  # Match the rounding from Grader.__call__
            except Exception as e:
                logger.error(f"Unexpected error during h5ad grading: {e}")
                raise e
        else:
            logger.warning(f"Unsupported file format: {file_extension}. Expected .csv or .h5ad.")
    else:
        logger.warning(
            f"Invalid submission file: {path_to_submission}. Please check that the file exists."
        )

    valid_submission = score is not None
    task_leaderboard = get_leaderboard(task)
    rank_info = task.grader.rank_score(score, task_leaderboard)
    is_lower_better = task.grader.is_lower_better(task_leaderboard)

    # Check human baselines if available
    beats_human = None
    human_percentile = None
    if score is not None:
        human_baselines_path = task.public_dir / "human_baselines.csv"
        if human_baselines_path.exists():
            try:
                human_df = pd.read_csv(human_baselines_path)
                if not human_df.empty and "score" in human_df.columns:
                    beats_human, human_percentile = calculate_human_performance_metrics(
                        score, human_df, task.grader.is_lower_better(task_leaderboard)
                    )
                    logger.debug(
                        f"Human baseline comparison: beats_human={beats_human}, percentile={human_percentile}"
                    )
            except Exception as e:
                logger.warning(f"Failed to load human baselines for task '{task.id}': {e}")

    return TaskReport(
        task_id=task.id,
        score=score,
        gold_threshold=rank_info["gold_threshold"],
        silver_threshold=rank_info["silver_threshold"],
        bronze_threshold=rank_info["bronze_threshold"],
        median_threshold=rank_info["median_threshold"],
        any_medal=rank_info["gold_medal"] or rank_info["silver_medal"] or rank_info["bronze_medal"],
        gold_medal=rank_info["gold_medal"],
        silver_medal=rank_info["silver_medal"],
        bronze_medal=rank_info["bronze_medal"],
        above_median=rank_info["above_median"],
        submission_exists=submission_exists,
        valid_submission=valid_submission,
        is_lower_better=is_lower_better,
        created_at=datetime.now(),
        submission_path=str(path_to_submission),
        leaderboard_position=rank_info["leaderboard_position"],
        beats_human=beats_human,
        human_percentile=human_percentile,
    )


def validate_submission(submission: Path, task: Task) -> tuple[bool, str]:
    """
    Validates a submission for the given task by actually running the task grader.
    This is designed for end users, not developers (we assume that the task grader is
    correctly implemented and use that for validating the submission, not the other way around).
    """
    if not submission.is_file():
        return False, f"Submission invalid! Submission file {submission} does not exist."

    # Support multiple file formats: CSV, h5ad, and potentially others
    supported_formats = [".csv", ".h5ad"]
    file_extension = submission.suffix.lower()
    if file_extension not in supported_formats:
        return (
            False,
            f"Submission invalid! Submission file must be one of {supported_formats}, got {file_extension}.",
        )

    if not is_dataset_prepared(task, grading_only=True):
        raise ValueError(
            f"Dataset for task `{task.id}` is not prepared! "
            f"Please run `biomlbench prepare -t {task.id}` to prepare the dataset."
        )

    try:
        # Use the same logic as grade_submission for different file formats
        if file_extension == ".csv":
            task.grader.grade_fn(read_csv(submission), read_csv(task.answers))
        elif file_extension == ".h5ad":
            # For h5ad files, pass paths directly to grade_fn
            task.grader.grade_fn(submission, Path(task.answers))
        else:
            # This shouldn't happen given our validation above, but handle it gracefully
            raise ValueError(f"Unsupported file format for validation: {file_extension}")
    except Exception as e:
        return (
            False,
            f"Submission invalid! The attempt to grade the submission has resulted in the following error message:\n{e}",
        )

    return True, "Submission is valid."


def aggregate_reports(task_reports: list[TaskReport]) -> dict:
    """
    Builds the summary report from a list of task reports.
    """

    total_gold_medals = sum(report.gold_medal for report in task_reports)
    total_silver_medals = sum(report.silver_medal for report in task_reports)
    total_bronze_medals = sum(report.bronze_medal for report in task_reports)
    total_above_median = sum(report.above_median for report in task_reports)
    total_submissions = sum(report.submission_exists for report in task_reports)
    total_valid_submissions = sum(report.valid_submission for report in task_reports)
    total_beats_human = sum(1 for report in task_reports if report.beats_human)

    # Calculate position statistics
    positions = [report.leaderboard_position for report in task_reports if report.leaderboard_position is not None]
    avg_position = sum(positions) / len(positions) if positions else None
    best_position = min(positions) if positions else None
    worst_position = max(positions) if positions else None

    summary_report = {
        "total_runs": int(len(task_reports)),
        "total_runs_with_submissions": int(total_submissions),
        "total_valid_submissions": int(total_valid_submissions),
        "total_medals": int(total_gold_medals + total_silver_medals + total_bronze_medals),
        "total_gold_medals": int(total_gold_medals),
        "total_silver_medals": int(total_silver_medals),
        "total_bronze_medals": int(total_bronze_medals),
        "total_above_median": int(total_above_median),
        "total_beats_human": int(total_beats_human),
        "avg_leaderboard_position": round(avg_position, 1) if avg_position is not None else None,
        "best_leaderboard_position": int(best_position) if best_position is not None else None,
        "worst_leaderboard_position": int(worst_position) if worst_position is not None else None,
        "task_reports": [tr.to_dict() for tr in task_reports],
    }

    return summary_report


def calculate_human_performance_metrics(
    agent_score: float, human_df: pd.DataFrame, is_lower_better: bool
) -> tuple[bool, float]:
    """
    Calculate how an agent's performance compares to human baselines.

    Args:
        agent_score: The agent's score to compare
        human_df: DataFrame with human baseline scores (must have 'score' column)
        is_lower_better: Whether lower scores are better for this metric

    Returns:
        Tuple of (beats_human, human_percentile) where:
        - beats_human: True if agent beats median human performance
        - human_percentile: Percentile of human performance that agent achieves (0-100)
    """
    human_scores = human_df["score"].astype(float)

    if human_scores.empty:
        return None, None

    # Calculate median human performance for beats_human
    median_human = human_scores.median()

    if is_lower_better:
        beats_human = agent_score < median_human
        # For lower-is-better: percentile = % of humans with score >= agent_score
        human_percentile = (human_scores >= agent_score).mean() * 100
    else:
        beats_human = agent_score > median_human
        # For higher-is-better: percentile = % of humans with score <= agent_score
        human_percentile = (human_scores <= agent_score).mean() * 100

    return beats_human, round(human_percentile, 1)
