import io
import shutil
import zipfile
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd
import requests

from biomlbench.data_sources.base import DataSource, DataSourceError
from biomlbench.data_sources.factory import register_data_source
from biomlbench.utils import get_logger

logger = get_logger(__name__)


@register_data_source("proteingym-dms")
class ProteinGymDMSDataSource(DataSource):
    """
    Data source for the ProteinGym deep mutational scanning (DMS) dataset.
    """

    DATASET_URL = "https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/cv_folds_singles_substitutions.zip"
    LEADERBOARD_URL = "https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/refs/heads/main/benchmarks/DMS_supervised/substitutions/Spearman/DMS_substitutions_Spearman_DMS_level.csv"
    METADATA_URL = "https://zenodo.org/records/15293562/files/DMS_substitutions.csv"

    # Class-level cache variables to share across all instances
    _cached_zip_path = None
    _cached_leaderboard = None
    _cached_metadata = None

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
        Download ProteinGym DMS substitutions data and extract to a shared location.
        The zip file is downloaded once and extracted once to a shared directory that all tasks can access.

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
        if self.__class__._cached_zip_path is None:
            shared_data_dir = data_dir.parent.parent
            zip_path = shared_data_dir / "cv_folds_singles_substitutions.zip"
            metadata_path = shared_data_dir / "DMS_substitutions.csv"
            extracted_dir = shared_data_dir / "cv_folds_singles_substitutions"
        else:
            # Subsequent times: use the cached shared location
            shared_data_dir = self.__class__._cached_zip_path.parent
            zip_path = self.__class__._cached_zip_path
            metadata_path = shared_data_dir / "DMS_substitutions.csv"
            extracted_dir = shared_data_dir / "cv_folds_singles_substitutions"

        # Download zip file if not cached
        if not (self.__class__._cached_zip_path and self.__class__._cached_zip_path.exists()):
            if not zip_path.exists():
                try:
                    # Ensure shared data directory exists
                    shared_data_dir.mkdir(parents=True, exist_ok=True)

                    logger.info(f"Downloading ProteinGym DMS data from {self.DATASET_URL}...")

                    # Download the zip file
                    response = requests.get(self.DATASET_URL, stream=True)
                    response.raise_for_status()

                    # Save the zip file
                    with open(zip_path, "wb") as f:
                        for chunk in response.iter_content(chunk_size=8192):
                            f.write(chunk)

                    logger.info(f"Successfully downloaded ProteinGym DMS data to {zip_path}")
                    self.__class__._cached_zip_path = zip_path

                except requests.RequestException as e:
                    raise DataSourceError(
                        f"Failed to download ProteinGym DMS data: {e}", source_type="proteingym"
                    ) from e
                except Exception as e:
                    raise DataSourceError(
                        f"Unexpected error downloading ProteinGym DMS data: {e}",
                        source_type="proteingym",
                    ) from e
            else:
                self.__class__._cached_zip_path = zip_path
        else:
            # Use cached zip file
            logger.info(f"Using cached ProteinGym zip: {self.__class__._cached_zip_path}")
            zip_path = self.__class__._cached_zip_path

        # Extract the zip file to shared location
        if not extracted_dir.exists():
            logger.info(f"Extracting ProteinGym data to shared location: {extracted_dir}...")
            try:
                with zipfile.ZipFile(zip_path, "r") as zip_ref:
                    zip_ref.extractall(shared_data_dir)
                logger.info(f"Successfully extracted ProteinGym data to {extracted_dir}")
            except Exception as e:
                raise DataSourceError(
                    f"Failed to extract ProteinGym data: {e}", source_type="proteingym"
                ) from e

        # Download metadata to shared location if not already present
        if not metadata_path.exists():
            try:
                logger.info(f"Downloading ProteinGym metadata to shared location...")
                response = requests.get(self.METADATA_URL)
                response.raise_for_status()

                with open(metadata_path, "w", encoding="utf-8") as f:
                    f.write(response.text)

                logger.info(f"Successfully downloaded metadata to {metadata_path}")
            except Exception as e:
                raise DataSourceError(
                    f"Failed to download ProteinGym metadata: {e}", source_type="proteingym"
                ) from e

        # Copy the specific dataset CSV file to the task's raw directory
        dataset_name = source_config["benchmark_id"]
        source_csv = extracted_dir / f"{dataset_name}.csv"
        target_csv = data_dir / f"{dataset_name}.csv"

        # Ensure the source CSV file exists
        if not source_csv.exists():
            raise DataSourceError(
                f"Dataset CSV file not found in extracted directory: {source_csv}",
                source_type="proteingym",
            )

        if not target_csv.exists():
            logger.info(f"Copying dataset CSV from shared location to task directory: {target_csv}")
            shutil.copy2(source_csv, target_csv)
        else:
            logger.info(f"Dataset CSV already exists in task directory: {target_csv}")

        # Return the task's raw directory (which now contains the specific CSV file)
        return data_dir

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

        # Check class-level cache first
        if self.__class__._cached_leaderboard is not None:
            logger.info("Using cached ProteinGym leaderboard")
            leaderboard_df = self.__class__._cached_leaderboard
        else:
            try:
                # Download and parse leaderboard CSV
                leaderboard_df = self._download_and_parse_leaderboard_csv()
                if leaderboard_df is not None and not leaderboard_df.empty:
                    self.__class__._cached_leaderboard = leaderboard_df
                else:
                    # If CSV parsing fails, return empty DataFrame
                    logger.warning("CSV parsing failed, returning empty leaderboard")
                    empty_leaderboard = pd.DataFrame(
                        columns=["teamName", "score", "submissionDate"]
                    )
                    self.__class__._cached_leaderboard = empty_leaderboard
                    return empty_leaderboard

            except Exception as e:
                logger.error(f"Error downloading/parsing ProteinGym leaderboard: {e}")
                empty_leaderboard = pd.DataFrame(columns=["teamName", "score", "submissionDate"])
                self.__class__._cached_leaderboard = empty_leaderboard
                return empty_leaderboard

        # If dataset_name is specified, filter to that specific dataset
        dataset_name = source_config["benchmark_id"]
        if dataset_name not in leaderboard_df.index:
            raise ValueError(f"Dataset {dataset_name} not found in leaderboard")

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
        logger.info(f"Filtered leaderboard to dataset: {dataset_name}")
        return dataset_df

    def _download_and_parse_leaderboard_csv(self) -> Optional[pd.DataFrame]:
        """
        Download and parse the ProteinGym leaderboard CSV file.

        Returns:
            DataFrame with leaderboard data or None if parsing fails
        """
        try:
            logger.info(f"Downloading ProteinGym leaderboard from {self.LEADERBOARD_URL}...")

            # Download the CSV file
            response = requests.get(self.LEADERBOARD_URL)
            response.raise_for_status()

            # Parse CSV content
            csv_content = response.text
            df = pd.read_csv(io.StringIO(csv_content))

            logger.info(f"Successfully downloaded CSV with shape: {df.shape}")
            logger.debug(f"CSV columns: {list(df.columns)}")
            logger.debug(f"First few rows:\n{df.head()}")

            # The CSV should have dataset names as the first column and team/model names as other columns
            # Set the first column as index (dataset names)
            if len(df.columns) > 0:
                df.set_index(df.columns[0], inplace=True)
                logger.info(
                    f"Parsed leaderboard with {len(df)} datasets and {len(df.columns)} teams/models"
                )
                return df
            else:
                logger.warning("CSV file appears to be empty or malformed")
                return None

        except requests.RequestException as e:
            logger.error(f"Failed to download CSV file: {e}")
            raise DataSourceError(
                f"Failed to download ProteinGym leaderboard CSV: {e}", source_type="proteingym"
            ) from e
        except pd.errors.EmptyDataError as e:
            logger.error(f"CSV file is empty: {e}")
            raise DataSourceError(
                f"ProteinGym leaderboard CSV is empty: {e}", source_type="proteingym"
            ) from e
        except pd.errors.ParserError as e:
            logger.error(f"Failed to parse CSV file: {e}")
            raise DataSourceError(
                f"Failed to parse ProteinGym leaderboard CSV: {e}", source_type="proteingym"
            ) from e
        except Exception as e:
            logger.error(f"Unexpected error downloading/parsing CSV: {e}")
            raise DataSourceError(
                f"Unexpected error with ProteinGym leaderboard CSV: {e}", source_type="proteingym"
            ) from e
