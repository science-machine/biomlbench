"""
Data preparation and management for BioML-bench.

This module provides high-level functions for downloading and preparing
datasets from various sources (Kaggle, Polaris, etc.) using a pluggable
data source architecture.
"""

import functools
import hashlib
import os
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Optional

import diskcache as dc
import pandas as pd
import yaml
from tqdm.auto import tqdm

from biomlbench.data_sources import DataSourceError, DataSourceFactory
from biomlbench.registry import Task
from biomlbench.utils import (
    extract,
    get_diff,
    get_logger,
    get_path_to_callable,
    is_empty,
    load_yaml,
)

logger = get_logger(__name__)
cache = dc.Cache("cache", size_limit=2**26)  # 64 MB


def create_prepared_dir(task: Task) -> None:
    """Create public and private directories for a task."""
    task.public_dir.mkdir(exist_ok=True, parents=True)
    task.private_dir.mkdir(exist_ok=True, parents=True)


def download_and_prepare_dataset(
    task: Task,
    keep_raw: bool = True,
    overwrite_checksums: bool = False,
    overwrite_leaderboard: bool = False,
    skip_verification: bool = False,
) -> None:
    """
    Download and prepare a dataset using the appropriate data source.

    Args:
        task: Task to prepare
        keep_raw: Whether to keep raw downloaded data
        overwrite_checksums: Whether to overwrite existing checksums
        overwrite_leaderboard: Whether to overwrite existing leaderboard
        skip_verification: Whether to skip checksum verification
    """

    assert is_valid_prepare_fn(
        task.prepare_fn
    ), f"Provided `prepare_fn` doesn't take arguments `raw`, `private` and `public`!"

    # Ensure leaderboard exists
    ensure_leaderboard_exists(task, force=overwrite_leaderboard)

    task_dir = task.raw_dir.parent
    task.raw_dir.mkdir(exist_ok=True, parents=True)
    create_prepared_dir(task)

    # Get data source configuration
    source_config = getattr(task, "data_source", None)
    if source_config is None:
        raise ValueError(f"No data_source configuration found for task '{task.id}'.")

    # Create appropriate data source
    source_type = source_config.get("type")
    try:
        data_source = DataSourceFactory.create(source_type)
    except Exception as e:
        raise ValueError(f"Failed to create data source '{source_type}': {e}") from e

    # Download data using the data source
    try:
        downloaded_path = data_source.download(source_config, task.raw_dir)

        # Handle zip file extraction for sources that provide zip files (like Kaggle)
        # Skip extraction for data sources that handle it internally (like ProteinGym)
        if downloaded_path and downloaded_path.suffix == ".zip":
            # Check if the zip file exists and hasn't been extracted yet
            # Look for the extracted directory that would be created by the zip
            zip_name = downloaded_path.stem  # Remove .zip extension
            extracted_dir = task.raw_dir / zip_name

            if not extracted_dir.exists():
                logger.info(f"Extracting `{downloaded_path}` to `{task.raw_dir}`...")
                extract(downloaded_path, task.raw_dir, recursive=False)
                logger.info(f"Extracted successfully.")
            else:
                logger.info(f"Zip file already extracted to `{extracted_dir}`")

    except DataSourceError as e:
        raise ValueError(f"Failed to download data for task '{task.id}': {e}") from e

    # Handle checksums if needed
    if overwrite_checksums or not skip_verification:
        actual_checksums = {}

        # Only include zip checksum if we have a zip file (and it's not handled internally)
        if downloaded_path and downloaded_path.suffix == ".zip":
            logger.info(f"Generating checksum for `{downloaded_path}`...")
            actual_zip_checksum = get_checksum(downloaded_path)
            actual_checksums["zip"] = actual_zip_checksum

            if task.checksums.is_file() and not overwrite_checksums:
                expected_checksums = load_yaml(task.checksums)
                expected_zip_checksum = expected_checksums.get("zip")

                if expected_zip_checksum and actual_zip_checksum != expected_zip_checksum:
                    raise ValueError(
                        f"Checksum for `{downloaded_path}` does not match the expected checksum! "
                        f"Expected `{expected_zip_checksum}` but got `{actual_zip_checksum}`."
                    )

                logger.info(f"Checksum for `{downloaded_path}` matches the expected checksum.")

    # Run task-specific preparation
    if not is_dataset_prepared(task) or overwrite_checksums:
        if task.public_dir.parent.exists() and overwrite_checksums:
            logger.info(
                f"Removing the existing prepared data directory for `{task.id}` since "
                "`overwrite_checksums` is set to `True`..."
            )
            shutil.rmtree(task.public_dir.parent)
            create_prepared_dir(task)

        logger.info(
            f"Preparing the dataset using `{task.prepare_fn.__name__}` from "
            f"`{get_path_to_callable(task.prepare_fn)}`..."
        )

        task.prepare_fn(
            raw=task.raw_dir,
            public=task.public_dir,
            private=task.private_dir,
        )

        logger.info(f"Data for task `{task.id}` prepared successfully.")

    # Save task description
    with open(task.public_dir / "description.md", "w") as f:
        f.write(task.description)

    # Prepare human baselines
    prepare_human_baselines(task, force=overwrite_checksums)

    # Generate final checksums
    if overwrite_checksums or not skip_verification:
        logger.info(f"Generating checksums for files in `{task_dir}`...")

        actual_checksums.update(
            {
                "public": generate_checksums(task.public_dir),
                "private": generate_checksums(task.private_dir),
            }
        )

        if not task.checksums.is_file() or overwrite_checksums:
            with open(task.checksums, "w") as file:
                yaml.dump(actual_checksums, file, default_flow_style=False)

            logger.info(f"Checksums for `{task.id}` saved to `{task.checksums}`.")

        if task.checksums.is_file() and not skip_verification:
            expected_checksums = load_yaml(task.checksums)

            if actual_checksums != expected_checksums:
                logger.error(f"Checksums do not match for `{task.id}`!")

                diff = get_diff(
                    actual_checksums,
                    expected_checksums,
                    fromfile="actual_checksums",
                    tofile="expected_checksums",
                )

                raise ValueError(f"Checksums do not match for `{task.id}`!\n{diff}")

            logger.info(f"Checksums for files in `{task_dir}` match the expected checksums.")

    # Clean up raw data if requested
    if not keep_raw and downloaded_path:
        logger.info(f"Removing the raw data directory for `{task.id}`...")
        if downloaded_path.is_file():
            downloaded_path.unlink()
        elif downloaded_path.is_dir():
            shutil.rmtree(downloaded_path)

    # Final validation
    assert task.public_dir.is_dir(), f"Public data directory doesn't exist."
    assert task.private_dir.is_dir(), f"Private data directory doesn't exist."
    assert not is_empty(task.public_dir), f"Public data directory is empty!"
    assert not is_empty(task.private_dir), f"Private data directory is empty!"


def is_dataset_prepared(task: Task, grading_only: bool = False) -> bool:
    """Checks if the task has non-empty `public` and `private` directories with the expected files."""

    assert isinstance(task, Task), f"Expected input to be of type `Task` but got {type(task)}."

    public = task.public_dir
    private = task.private_dir

    if not grading_only:
        if not public.is_dir():
            logger.warning("Public directory does not exist.")
            return False
        if is_empty(public):
            logger.warning("Public directory is empty.")
            return False

    if not private.is_dir():
        logger.warning("Private directory does not exist.")
        return False
    if is_empty(private):
        logger.warning("Private directory is empty.")
        return False

    if not task.answers.is_file():
        logger.warning("Answers file does not exist.")
        return False

    if not task.sample_submission.is_file() and not grading_only:
        logger.warning("Sample submission file does not exist.")
        return False

    # Check for leaderboard - critical for grading and medal calculations
    if not task.leaderboard.is_file():
        logger.warning("Leaderboard file does not exist.")
        return False

    return True


def ensure_leaderboard_exists(task: Task, force: bool = False) -> Path:
    """
    Ensures the leaderboard for a given task exists.

    Args:
        task: Task to ensure leaderboard for
        force: Whether to force download/update of leaderboard

    Returns:
        Path to the leaderboard file

    Raises:
        FileNotFoundError: If leaderboard cannot be found or created
    """
    leaderboard_path = task.leaderboard

    if not force and leaderboard_path.exists():
        return leaderboard_path

    # Try to get leaderboard from data source
    source_config = getattr(task, "data_source", None)
    if source_config is None:
        # Fallback: assume Kaggle for backward compatibility
        source_config = {"type": "kaggle", "competition_id": task.id}

    source_type = source_config.get("type", "kaggle")

    try:
        data_source = DataSourceFactory.create(source_type)
        leaderboard_df = data_source.get_leaderboard(source_config)

        # Save leaderboard
        leaderboard_path.parent.mkdir(parents=True, exist_ok=True)
        leaderboard_df.to_csv(leaderboard_path, index=False)

        logger.info(f"Downloaded leaderboard for task `{task.id}` using {source_type} data source")
        return leaderboard_path

    except Exception as e:
        if not force and leaderboard_path.exists():
            # Fallback to existing file if download fails
            logger.warning(f"Failed to update leaderboard for task `{task.id}`: {e}")
            return leaderboard_path
        else:
            raise FileNotFoundError(
                f"Leaderboard not found locally for task `{task.id}` and could not be downloaded: {e}"
            ) from e


def is_valid_prepare_fn(preparer_fn: Any) -> bool:
    """Checks if the `preparer_fn` takes three arguments: `raw`, `public` and `private`, in that order."""

    import inspect

    try:
        sig = inspect.signature(preparer_fn)
    except (TypeError, ValueError):
        return False

    actual_params = list(sig.parameters.keys())
    expected_params = ["raw", "public", "private"]

    return actual_params == expected_params


def generate_checksums(
    target_dir: Path,
    exts: Optional[list[str]] = None,
    exclude: Optional[list[Path]] = None,
) -> dict:
    """
    Generate checksums for the files directly under the target directory with the specified extensions.

    Args:
        target_dir: directory to generate checksums for.
        exts: List of file extensions to generate checksums for.
        exclude: List of file paths to exclude from checksum generation.

    Returns:
        A dictionary of form file: checksum.
    """

    if exts is None:
        exts = ["csv", "json", "jsonl", "parquet", "bson", "h5ad"]

    if exclude is None:
        exclude = []

    checksums = {}

    for ext in exts:
        fpaths = target_dir.glob(f"*.{ext}")

        for fpath in fpaths:
            if not fpath.is_file():
                continue  # skip dirs named like `my/dir.csv/`

            if fpath in exclude:
                continue

            checksums[fpath.name] = get_checksum(fpath)

    return checksums


def get_last_modified(fpath: Path) -> datetime:
    """Return the last modified time of a file."""

    return datetime.fromtimestamp(fpath.stat().st_mtime)


def file_cache(fn: Callable) -> Callable:
    """A decorator that caches results of a function with a Path argument, invalidating the cache when the file is modified."""

    import inspect

    sig = inspect.signature(fn)
    params = list(sig.parameters.values())

    if not (len(params) == 1 and params[0].annotation is Path):
        raise NotImplementedError("Only functions with a single `Path` argument are supported.")

    # Use `functools.wraps` to preserve the function's metadata, like the name and docstring.
    # Query the cache, but with an additional `last_modified` argument in the key, which has the
    # side effect of invalidating the cache when the file is modified.
    @functools.wraps(fn)
    def wrapper(fpath: Path) -> Any:
        last_modified = get_last_modified(fpath)
        key = (fn.__name__, str(fpath), last_modified)

        if key not in cache:
            cache[key] = fn(fpath)

        return cache[key]

    return wrapper


@file_cache
def get_checksum(fpath: Path) -> str:
    """Compute MD5 checksum of a file."""

    assert fpath.is_file(), f"Expected a file at `{fpath}`, but it doesn't exist."

    hash_md5 = hashlib.md5()
    file_size = os.path.getsize(fpath)

    # only show progress bar for large files (> ~5 GB)
    show_progress = file_size > 5_000_000_000

    with open(fpath, "rb") as f:
        for chunk in tqdm(
            iter(lambda: f.read(4_096), b""),
            total=file_size // 4096,
            unit="B",
            unit_scale=True,
            disable=not show_progress,
        ):
            hash_md5.update(chunk)

    return hash_md5.hexdigest()


def get_leaderboard(task: Task) -> pd.DataFrame:
    """Load leaderboard data for a task."""
    leaderboard_path = task.leaderboard
    assert leaderboard_path.exists(), f"Leaderboard not found locally for task `{task.id}`."
    leaderboard_df = pd.read_csv(leaderboard_path)
    return leaderboard_df


def prepare_human_baselines(task: Task, force: bool = False) -> Optional[Path]:
    """
    Prepare human baseline data for a task.

    Args:
        task: Task to prepare human baselines for
        force: Whether to force re-download of human baselines

    Returns:
        Path to human baselines CSV file, or None if not available
    """
    human_baselines_path = task.public_dir / "human_baselines.csv"

    # Skip if already exists and not forcing
    if human_baselines_path.exists() and not force:
        logger.info(f"Human baselines already exist for task '{task.id}'")
        return human_baselines_path

    # Get data source configuration
    source_config = getattr(task, "data_source", None)
    if source_config is None:
        logger.debug(f"No data source configuration for task '{task.id}', skipping human baselines")
        return None

    # Create appropriate data source
    source_type = source_config.get("type")
    try:
        data_source = DataSourceFactory.create(source_type)
    except Exception as e:
        logger.warning(f"Failed to create data source '{source_type}' for human baselines: {e}")
        return None

    # Check if data source supports human baselines
    if not data_source.supports_human_baselines():
        logger.debug(f"Data source '{source_type}' does not support human baselines")
        return None

    # Extract human baselines
    try:
        human_baselines_df = data_source.get_human_baselines(source_config)

        if human_baselines_df is None or human_baselines_df.empty:
            logger.info(f"No human baselines available for task '{task.id}'")
            return None

        # Save human baselines
        human_baselines_df.to_csv(human_baselines_path, index=False)

        logger.info(f"Saved {len(human_baselines_df)} human baseline entries for task '{task.id}'")
        return human_baselines_path

    except Exception as e:
        logger.error(f"Failed to extract human baselines for task '{task.id}': {e}")
        raise e
