import io
import shutil
import zipfile
from pathlib import Path
from typing import Any, Dict, Literal, Optional

import pandas as pd
import requests

from biomlbench.data_sources.base import DataSource, DataSourceError
from biomlbench.data_sources.factory import register_data_source
from biomlbench.utils import get_logger

logger = get_logger(__name__)


@register_data_source("proteingym-dms")
class ProteinGymDMSDataSource(DataSource):
    """
    Data source for the ProteinGym deep mutational scanning (DMS) benchmark.
    """

    SINGLE_SUBSTITUTIONS_URL = "https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/cv_folds_singles_substitutions.zip"
    MULTI_SUBSTITUTIONS_URL = "https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/cv_folds_multiples_substitutions.zip"
    INDELS_URL = "https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/cv_folds_indels.zip"

    SUBSTITUTIONS_METADATA_URL = "https://zenodo.org/records/15293562/files/DMS_substitutions.csv"
    INDELS_METADATA_URL = "https://zenodo.org/records/15293562/files/DMS_indels.csv"

    SUBSTITUTIONS_LEADERBOARD_URL = "https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/refs/heads/main/benchmarks/DMS_supervised/substitutions/Spearman/DMS_substitutions_Spearman_DMS_level.csv"
    INDELS_LEADERBOARD_URL = "https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/refs/heads/main/benchmarks/DMS_supervised/indels/Spearman/DMS_indels_Spearman_DMS_level.csv"

    # Class-level cache variables to share across all instances
    _cached_single_substitutions_zip_path = None
    _cached_multi_substitutions_zip_path = None
    _cached_indels_zip_path = None
    _cached_substitutions_metadata_path = None
    _cached_indels_metadata_path = None
    _cached_substitutions_leaderboard = None
    _cached_indels_leaderboard = None

    def validate_config(self, source_config: Dict[str, Any]) -> bool:
        """
        Validate ProteinGym source configuration.

        Args:
            source_config: May contain optional 'dataset_name' key

        Returns:
            True if valid

        Raises:
            ValueError: If configuration is invalid
        """
        # Optional dataset_name for indexing specific dataset in leaderboard
        if "benchmark_id" not in source_config:
            raise ValueError("ProteinGym 'benchmark_id' is required")

        dataset_name = source_config["benchmark_id"]
        if not isinstance(dataset_name, str) or not dataset_name.strip():
            raise ValueError("ProteinGym 'benchmark_id' must be a non-empty string")

        return True

    def download(self, source_config: Dict[str, Any], data_dir: Path) -> Path | None:
        """
        Download ProteinGym DMS data (single substitutions, multi substitutions, and indels)
        and extract to a shared location. The zip files are downloaded once and extracted
        once to a shared directory that all tasks can access.

        Args:
            source_config: Configuration for the data source, must contain 'benchmark_id'
            data_dir: Directory to download data to (not used for extraction)

        Returns:
            Path to the task's raw directory (containing the specific dataset CSV file).

        Raises:
            DataSourceError: If download fails
        """
        # Validate config first
        self.validate_config(source_config)

        # Define shared paths (use the proteingym-dms level as the shared location)
        shared_data_dir = data_dir.parent.parent

        # Define zip file paths
        single_substitutions_zip_path = shared_data_dir / "cv_folds_singles_substitutions.zip"
        multi_substitutions_zip_path = shared_data_dir / "cv_folds_multiples_substitutions.zip"
        indels_zip_path = shared_data_dir / "cv_folds_indels.zip"

        # Define metadata file paths
        substitutions_metadata_path = shared_data_dir / "DMS_substitutions.csv"
        indels_metadata_path = shared_data_dir / "DMS_indels.csv"

        # Define extracted directory paths
        single_substitutions_dir = shared_data_dir / "cv_folds_singles_substitutions"
        multi_substitutions_dir = shared_data_dir / "cv_folds_multiples_substitutions"
        indels_dir = shared_data_dir / "cv_folds_indels"

        # Download all three datasets if not already cached
        self._download_dataset_if_needed(
            self.SINGLE_SUBSTITUTIONS_URL, single_substitutions_zip_path, "single"
        )
        self._download_dataset_if_needed(
            self.MULTI_SUBSTITUTIONS_URL, multi_substitutions_zip_path, "multi"
        )
        self._download_dataset_if_needed(self.INDELS_URL, indels_zip_path, "indels")

        # Extract all datasets if not already extracted
        self._extract_dataset_if_needed(
            single_substitutions_zip_path, single_substitutions_dir, "single substitutions"
        )
        self._extract_dataset_if_needed(
            multi_substitutions_zip_path, multi_substitutions_dir, "multi substitutions"
        )
        self._extract_dataset_if_needed(indels_zip_path, indels_dir, "indels")

        # Download metadata files to shared location if not already present
        self._download_metadata_if_needed(
            self.SUBSTITUTIONS_METADATA_URL, substitutions_metadata_path, "substitutions"
        )
        self._download_metadata_if_needed(self.INDELS_METADATA_URL, indels_metadata_path, "indels")

        # Find the specific dataset CSV file in one of the extracted directories
        dataset_name = source_config["benchmark_id"]
        source_csv = None

        # Check in single substitutions first
        potential_source_csv = single_substitutions_dir / f"{dataset_name}.csv"
        if potential_source_csv.exists():
            source_csv = potential_source_csv
            logger.info(f"Found dataset {dataset_name} in single substitutions directory")
        else:
            # Check in multi substitutions
            potential_source_csv = multi_substitutions_dir / f"{dataset_name}.csv"
            if potential_source_csv.exists():
                source_csv = potential_source_csv
                logger.info(f"Found dataset {dataset_name} in multi substitutions directory")
            else:
                # Check in indels
                potential_source_csv = indels_dir / f"{dataset_name}.csv"
                if potential_source_csv.exists():
                    source_csv = potential_source_csv
                    logger.info(f"Found dataset {dataset_name} in indels directory")

        if source_csv is None:
            raise DataSourceError(
                f"Dataset CSV file '{dataset_name}.csv' not found in any extracted directory. "
                f"Checked: {single_substitutions_dir}, {multi_substitutions_dir}, {indels_dir}",
                source_type="proteingym",
            )

        # Copy the specific dataset CSV file to the task's raw directory
        target_csv = data_dir / f"{dataset_name}.csv"

        if not target_csv.exists():
            logger.info(f"Copying dataset CSV from shared location to task directory: {target_csv}")
            shutil.copy2(source_csv, target_csv)
        else:
            logger.info(f"Dataset CSV already exists in task directory: {target_csv}")

        # Return the task's raw directory (which now contains the specific CSV file)
        return data_dir

    def _download_dataset_if_needed(
        self, url: str, zip_path: Path, dataset_type: Literal["single", "multi", "indels"]
    ) -> None:
        """
        Download a dataset zip file if it doesn't already exist.

        Args:
            url: URL to download from
            zip_path: Path where to save the zip file
            dataset_type: Human-readable name for logging

        Raises:
            DataSourceError: If download fails
        """
        if not zip_path.exists():
            try:
                # Ensure shared data directory exists
                zip_path.parent.mkdir(parents=True, exist_ok=True)

                logger.info(f"Downloading ProteinGym {dataset_type} data from {url}...")

                # Download the zip file
                response = requests.get(url, stream=True)
                response.raise_for_status()

                # Save the zip file
                with open(zip_path, "wb") as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)

                logger.info(f"Successfully downloaded ProteinGym {dataset_type} data to {zip_path}")

                # Update class-level cache
                if dataset_type == "single":
                    self.__class__._cached_single_substitutions_zip_path = zip_path
                elif dataset_type == "multi":
                    self.__class__._cached_multi_substitutions_zip_path = zip_path
                elif dataset_type == "indels":
                    self.__class__._cached_indels_zip_path = zip_path

            except requests.RequestException as e:
                raise DataSourceError(
                    f"Failed to download ProteinGym {dataset_type} data: {e}",
                    source_type="proteingym",
                ) from e
            except Exception as e:
                raise DataSourceError(
                    f"Unexpected error downloading ProteinGym {dataset_type} data: {e}",
                    source_type="proteingym",
                ) from e
        else:
            logger.info(f"Using cached ProteinGym {dataset_type} zip: {zip_path}")

            # Update class-level cache even if using existing file
            if dataset_type == "single":
                self.__class__._cached_single_substitutions_zip_path = zip_path
            elif dataset_type == "multi":
                self.__class__._cached_multi_substitutions_zip_path = zip_path
            elif dataset_type == "indels":
                self.__class__._cached_indels_zip_path = zip_path

    def _extract_dataset_if_needed(
        self, zip_path: Path, extracted_dir: Path, dataset_type: str
    ) -> None:
        """
        Extract a dataset zip file if it hasn't already been extracted.

        Args:
            zip_path: Path to the zip file
            extracted_dir: Directory to extract to
            dataset_type: Human-readable name for logging

        Raises:
            DataSourceError: If extraction fails
        """
        if not extracted_dir.exists():
            logger.info(
                f"Extracting ProteinGym {dataset_type} data to shared location: {extracted_dir}..."
            )
            try:
                with zipfile.ZipFile(zip_path, "r") as zip_ref:
                    zip_ref.extractall(extracted_dir.parent)
                logger.info(
                    f"Successfully extracted ProteinGym {dataset_type} data to {extracted_dir}"
                )
            except Exception as e:
                raise DataSourceError(
                    f"Failed to extract ProteinGym {dataset_type} data: {e}",
                    source_type="proteingym",
                ) from e
        else:
            logger.info(f"Using cached extracted ProteinGym {dataset_type} data: {extracted_dir}")

    def _download_metadata_if_needed(
        self, url: str, metadata_path: Path, metadata_type: Literal["substitutions", "indels"]
    ) -> None:
        """
        Download a metadata file if it doesn't already exist.

        Args:
            url: URL to download from
            metadata_path: Path where to save the metadata file
            metadata_type: Human-readable name for logging

        Raises:
            DataSourceError: If download fails
        """
        if not metadata_path.exists():
            try:
                logger.info(f"Downloading ProteinGym {metadata_type} to shared location...")
                response = requests.get(url)
                response.raise_for_status()

                with open(metadata_path, "w", encoding="utf-8") as f:
                    f.write(response.text)

                logger.info(f"Successfully downloaded {metadata_type} to {metadata_path}")

                # Update class-level cache
                if metadata_type == "substitutions":
                    self.__class__._cached_substitutions_metadata_path = metadata_path
                elif metadata_type == "indels":
                    self.__class__._cached_indels_metadata_path = metadata_path

            except Exception as e:
                raise DataSourceError(
                    f"Failed to download ProteinGym {metadata_type}: {e}", source_type="proteingym"
                ) from e
        else:
            logger.info(f"Using cached ProteinGym {metadata_type}: {metadata_path}")

            # Update class-level cache even if using existing file
            if metadata_type == "substitutions":
                self.__class__._cached_substitutions_metadata_path = metadata_path
            elif metadata_type == "indels":
                self.__class__._cached_indels_metadata_path = metadata_path

    def _determine_leaderboard_type(self, dataset_name: str) -> Literal["substitutions", "indels"]:
        """
        Determine which leaderboard to use based on the dataset name.

        Args:
            dataset_name: Name of the dataset

        Returns:
            Leaderboard type ("substitutions" or "indels")

        Raises:
            ValueError: If dataset is not found in either metadata file
        """
        # Check if the dataset name contains "indels"
        if "indels" in dataset_name.lower():
            return "indels"
        else:
            return "substitutions"

    def get_leaderboard(self, source_config: Dict[str, Any]) -> pd.DataFrame:
        """
        Get leaderboard from ProteinGym.

        Args:
            source_config: Configuration for the data source, must contain 'benchmark_id'

        Returns:
            DataFrame with columns: teamName, score, submissionDate
            If dataset_name is specified, returns only that dataset's row

        Raises:
            DataSourceError: If leaderboard cannot be retrieved
        """
        # Validate config first
        self.validate_config(source_config)

        # Determine which leaderboard to use based on dataset type
        dataset_name = source_config["benchmark_id"]
        leaderboard_type = self._determine_leaderboard_type(dataset_name)

        # Check class-level cache first
        cache_var = f"_cached_{leaderboard_type}_leaderboard"
        if getattr(self.__class__, cache_var) is not None:
            logger.info(f"Using cached ProteinGym {leaderboard_type} leaderboard")
            leaderboard_df = getattr(self.__class__, cache_var)
        else:
            try:
                # Download and parse appropriate leaderboard CSV
                leaderboard_df = self._download_and_parse_leaderboard_csv(leaderboard_type)
                if leaderboard_df is not None and not leaderboard_df.empty:
                    setattr(self.__class__, cache_var, leaderboard_df)
                else:
                    # If CSV parsing fails, return empty DataFrame
                    logger.warning(
                        f"CSV parsing failed for {leaderboard_type} leaderboard, returning empty leaderboard"
                    )
                    empty_leaderboard = pd.DataFrame(
                        columns=["teamName", "score", "submissionDate"]
                    )
                    setattr(self.__class__, cache_var, empty_leaderboard)
                    return empty_leaderboard

            except Exception as e:
                logger.error(
                    f"Error downloading/parsing ProteinGym {leaderboard_type} leaderboard: {e}"
                )
                empty_leaderboard = pd.DataFrame(columns=["teamName", "score", "submissionDate"])
                setattr(self.__class__, cache_var, empty_leaderboard)
                return empty_leaderboard

        # If dataset_name is specified, filter to that specific dataset
        if dataset_name not in leaderboard_df.index:
            raise ValueError(f"Dataset {dataset_name} not found in {leaderboard_type} leaderboard")

        # Extract the specific dataset row and convert to standard format
        dataset_row = leaderboard_df.loc[dataset_name]
        dataset_df = pd.DataFrame(
            {
                "teamName": dataset_row.index,
                "score": dataset_row.values,
                "submissionDate": "2024-01-01",  # Placeholder date
            }
        )
        dataset_df.sort_values(by="score", ascending=False, inplace=True)
        logger.info(f"Filtered {leaderboard_type} leaderboard to dataset: {dataset_name}")
        return dataset_df

    def _download_and_parse_leaderboard_csv(
        self, leaderboard_type: Literal["substitutions", "indels"]
    ) -> Optional[pd.DataFrame]:
        """
        Download and parse the ProteinGym leaderboard CSV file.

        Args:
            leaderboard_type: Type of leaderboard to download ("substitutions" or "indels")

        Returns:
            DataFrame with leaderboard data or None if parsing fails
        """
        # Determine the correct URL based on leaderboard type
        if leaderboard_type == "substitutions":
            url = self.SUBSTITUTIONS_LEADERBOARD_URL
        elif leaderboard_type == "indels":
            url = self.INDELS_LEADERBOARD_URL
        else:
            raise ValueError(f"Unknown leaderboard type: {leaderboard_type}")

        try:
            logger.info(f"Downloading ProteinGym {leaderboard_type} leaderboard from {url}...")

            # Download the CSV file
            response = requests.get(url)
            response.raise_for_status()

            # Parse CSV content
            csv_content = response.text
            df = pd.read_csv(io.StringIO(csv_content))

            logger.info(f"Successfully downloaded {leaderboard_type} CSV with shape: {df.shape}")
            logger.debug(f"CSV columns: {list(df.columns)}")
            logger.debug(f"First few rows:\n{df.head()}")

            # The CSV should have dataset names as the first column and team/model names as other columns
            # Set the first column as index (dataset names)
            if len(df.columns) > 0:
                df.set_index(df.columns[0], inplace=True)
                logger.info(
                    f"Parsed {leaderboard_type} leaderboard with {len(df)} datasets and {len(df.columns)} teams/models"
                )
                return df
            else:
                logger.warning(f"{leaderboard_type} CSV file appears to be empty or malformed")
                return None

        except requests.RequestException as e:
            logger.error(f"Failed to download {leaderboard_type} CSV file: {e}")
            raise DataSourceError(
                f"Failed to download ProteinGym {leaderboard_type} leaderboard CSV: {e}",
                source_type="proteingym",
            ) from e
        except pd.errors.EmptyDataError as e:
            logger.error(f"{leaderboard_type} CSV file is empty: {e}")
            raise DataSourceError(
                f"ProteinGym {leaderboard_type} leaderboard CSV is empty: {e}",
                source_type="proteingym",
            ) from e
        except pd.errors.ParserError as e:
            logger.error(f"Failed to parse {leaderboard_type} CSV file: {e}")
            raise DataSourceError(
                f"Failed to parse ProteinGym {leaderboard_type} leaderboard CSV: {e}",
                source_type="proteingym",
            ) from e
        except Exception as e:
            logger.error(f"Unexpected error downloading/parsing {leaderboard_type} CSV: {e}")
            raise DataSourceError(
                f"Unexpected error with ProteinGym {leaderboard_type} leaderboard CSV: {e}",
                source_type="proteingym",
            ) from e
