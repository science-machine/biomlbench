import argparse
import json
from pathlib import Path

from biomlbench.data import download_and_prepare_dataset, ensure_leaderboard_exists, prepare_human_baselines
from biomlbench.grade import grade_csv, grade_jsonl
from biomlbench.registry import registry
from biomlbench.utils import get_logger
from biomlbench.baselines import run_baseline

logger = get_logger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Runs agents on biomedical ML tasks.")
    subparsers = parser.add_subparsers(dest="command", help="Sub-command to run.")

    # Prepare sub-parser
    parser_prepare = subparsers.add_parser(
        name="prepare",
        help="Download and prepare tasks for the BioML-bench dataset.",
    )
    parser_prepare.add_argument(
        "-t",
        "--task-id",
        help=f"ID of the task to prepare. Valid options: {registry.list_task_ids()}",
        type=str,
        required=False,
    )
    parser_prepare.add_argument(
        "-a",
        "--all",
        help="Prepare all tasks.",
        action="store_true",
    )
    parser_prepare.add_argument(
        "--lite",
        help="Prepare all the low difficulty tasks (BioML-bench Lite).",
        action="store_true",
        required=False,
    )
    parser_prepare.add_argument(
        "-l",
        "--list",
        help="Prepare a list of tasks specified line by line in a text file.",
        type=str,
        required=False,
    )
    parser_prepare.add_argument(
        "--domain",
        help="Prepare all tasks for a specific biomedical domain (e.g., oncology, drug_discovery).",
        type=str,
        required=False,
    )
    parser_prepare.add_argument(
        "--task-type",
        help="Prepare all tasks of a specific type (e.g., medical_imaging, protein_engineering).",
        type=str,
        required=False,
    )
    parser_prepare.add_argument(
        "--keep-raw",
        help="Keep the raw task files after the task has been prepared.",
        action="store_true",
        required=False,
        default=False,
    )
    parser_prepare.add_argument(
        "--data-dir",
        help="Path to the directory where the data will be stored.",
        required=False,
        default=registry.get_data_dir(),
    )
    parser_prepare.add_argument(
        "--overwrite-checksums",
        help="[For Developers] Overwrite the checksums file for the task.",
        action="store_true",
        required=False,
        default=False,
    )
    parser_prepare.add_argument(
        "--overwrite-leaderboard",
        help="[For Developers] Overwrite the leaderboard file for the task.",
        action="store_true",
        required=False,
        default=False,
    )
    parser_prepare.add_argument(
        "--skip-verification",
        help="[For Developers] Skip the verification of the checksums.",
        action="store_true",
        required=False,
        default=False,
    )

    # Grade eval sub-parser
    parser_grade_eval = subparsers.add_parser(
        "grade",
        help="Grade a submission to the eval, comprising of several task submissions",
    )
    parser_grade_eval.add_argument(
        "--submission",
        help="Path to the JSONL file of submissions. Refer to README.md#submission-format for the required format.",
        type=str,
        required=True,
    )
    parser_grade_eval.add_argument(
        "--output-dir",
        help="Path to the directory where the evaluation metrics will be saved.",
        type=str,
        required=True,
    )
    parser_grade_eval.add_argument(
        "--data-dir",
        help="Path to the directory where the data used for grading is stored.",
        required=False,
        default=registry.get_data_dir(),
    )

    # Grade sample sub-parser
    parser_grade_sample = subparsers.add_parser(
        name="grade-sample",
        help="Grade a single sample (task) in the eval",
    )
    parser_grade_sample.add_argument(
        "submission",
        help="Path to the submission CSV file.",
        type=str,
    )
    parser_grade_sample.add_argument(
        "task_id",
        help=f"ID of the task to grade. Valid options: {registry.list_task_ids()}",
        type=str,
    )
    parser_grade_sample.add_argument(
        "--data-dir",
        help="Path to the directory where the data will be stored.",
        required=False,
        default=registry.get_data_dir(),
    )

    # Dev tools sub-parser
    parser_dev = subparsers.add_parser("dev", help="Developer tools for extending BioML-bench.")
    dev_subparsers = parser_dev.add_subparsers(dest="dev_command", help="Developer command to run.")

    # Set up 'download-leaderboard' under 'dev'
    parser_download_leaderboard = dev_subparsers.add_parser(
        "download-leaderboard",
        help="Download the leaderboard for a task.",
    )
    parser_download_leaderboard.add_argument(
        "-t",
        "--task-id",
        help=f"Name of the task to download the leaderboard for. Valid options: {registry.list_task_ids()}",
        type=str,
        required=False,
    )
    parser_download_leaderboard.add_argument(
        "--all",
        help="Download the leaderboard for all tasks.",
        action="store_true",
    )
    parser_download_leaderboard.add_argument(
        "--force",
        help="Force download the leaderboard, even if it already exists.",
        action="store_true",
    )

    # Set up 'prepare-human-baselines' under 'dev'
    parser_human_baselines = dev_subparsers.add_parser(
        "prepare-human-baselines",
        help="Prepare human baseline data for tasks.",
    )
    parser_human_baselines.add_argument(
        "-t",
        "--task-id",
        help=f"Name of the task to prepare human baselines for. Valid options: {registry.list_task_ids()}",
        type=str,
        required=False,
    )
    parser_human_baselines.add_argument(
        "--all",
        help="Prepare human baselines for all tasks.",
        action="store_true",
    )
    parser_human_baselines.add_argument(
        "--force",
        help="Force re-extraction of human baselines, even if they already exist.",
        action="store_true",
    )

    # Run baseline sub-parser
    parser_baseline = subparsers.add_parser(
        name="run-baseline",
        help="Run a baseline agent on a biomedical task",
    )
    parser_baseline.add_argument(
        "task_id",
        help=f"ID of the task to run baseline on. Valid options: {registry.list_task_ids()}",
        type=str,
    )
    parser_baseline.add_argument(
        "--baseline",
        help="Baseline to run (simple, random, or task-specific types like linear, rf, fingerprint). Use 'all' to run all available baselines for the task.",
        type=str,
        default="simple",
    )
    parser_baseline.add_argument(
        "--output-dir",
        help="Directory to save baseline submissions",
        type=str,
        default="baseline_submissions",
    )
    parser_baseline.add_argument(
        "--data-dir",
        help="Path to the directory where the data is stored.",
        required=False,
        default=registry.get_data_dir(),
    )
    parser_baseline.add_argument(
        "--seed",
        help="Random seed for reproducible baselines",
        type=int,
        default=42,
    )

    args = parser.parse_args()

    if args.command == "prepare":
        new_registry = registry.set_data_dir(Path(args.data_dir))

        if args.lite:
            tasks = [
                new_registry.get_task(task_id)
                for task_id in new_registry.get_lite_task_ids()
            ]
        elif args.all:
            tasks = [
                new_registry.get_task(task_id)
                for task_id in registry.list_task_ids()
            ]
        elif args.list:
            with open(args.list, "r") as f:
                task_ids = f.read().splitlines()
            tasks = [
                new_registry.get_task(task_id) for task_id in task_ids
            ]
        elif args.domain:
            task_ids = new_registry.get_tasks_by_domain(args.domain)
            tasks = [new_registry.get_task(task_id) for task_id in task_ids]
        elif args.task_type:
            task_ids = new_registry.get_tasks_by_type(args.task_type)
            tasks = [new_registry.get_task(task_id) for task_id in task_ids]
        else:
            if not args.task_id:
                parser_prepare.error(
                    "One of --lite, --all, --list, --domain, --task-type, or --task-id must be specified."
                )
            tasks = [new_registry.get_task(args.task_id)]

        for task in tasks:
            download_and_prepare_dataset(
                competition=task,  # Keep compatibility for now
                keep_raw=args.keep_raw,
                overwrite_checksums=args.overwrite_checksums,
                overwrite_leaderboard=args.overwrite_leaderboard,
                skip_verification=args.skip_verification,
            )
    
    if args.command == "grade":
        new_registry = registry.set_data_dir(Path(args.data_dir))
        submission = Path(args.submission)
        output_dir = Path(args.output_dir)
        grade_jsonl(submission, output_dir, new_registry)
    
    if args.command == "grade-sample":
        new_registry = registry.set_data_dir(Path(args.data_dir))
        task = new_registry.get_task(args.task_id)
        submission = Path(args.submission)
        report = grade_csv(submission, task)
        logger.info("Task report:")
        logger.info(json.dumps(report.to_dict(), indent=4))
    
    if args.command == "dev":
        if args.dev_command == "download-leaderboard":
            if args.all:
                for task_id in registry.list_task_ids():
                    task = registry.get_task(task_id)
                    ensure_leaderboard_exists(task, force=args.force)
            elif args.task_id:
                task = registry.get_task(args.task_id)
                ensure_leaderboard_exists(task, force=args.force)
            else:
                parser_download_leaderboard.error(
                    "Either --all or --task-id must be specified."
                )
        
        elif args.dev_command == "prepare-human-baselines":
            if args.all:
                for task_id in registry.list_task_ids():
                    task = registry.get_task(task_id)
                    prepare_human_baselines(task, force=args.force)
            elif args.task_id:
                task = registry.get_task(args.task_id)
                prepare_human_baselines(task, force=args.force)
            else:
                parser_human_baselines.error(
                    "Either --all or --task-id must be specified."
                )

    if args.command == "run-baseline":
        new_registry = registry.set_data_dir(Path(args.data_dir))
        task = new_registry.get_task(args.task_id)
        output_dir = Path(args.output_dir)
        
        run_baseline(
            task=task,
            baseline_type=args.baseline,
            output_dir=output_dir,
            seed=args.seed
        )


if __name__ == "__main__":
    main()
