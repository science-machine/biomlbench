import io
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd
import requests

from biomlbench.data_sources.base import DataSource, DataSourceError
from biomlbench.data_sources.factory import register_data_source
from biomlbench.utils import get_logger

logger = get_logger(__name__)


@register_data_source("proteingym")
class ProteinGymDMSDataSource(DataSource):
    """
    Data source for the ProteinGym deep mutational scanning (DMS) dataset.
    """

    DATASET_URL = "https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/cv_folds_singles_substitutions.zip"
    LEADERBOARD_URL = "https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/refs/heads/main/benchmarks/DMS_supervised/substitutions/Spearman/DMS_substitutions_Spearman_DMS_level.csv"
    METADATA_URL = "https://zenodo.org/records/15293562/files/DMS_substitutions.csv"

    def __init__(self):
        """Initialises the ProteinGym DMS data source without cached paths."""
        self._cached_zip_path = None
        self._cached_leaderboard = None
        self._cached_metadata = None

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
        Download ProteinGym DMS substitutions data as a zip file of all DMS datasets
        and a metadata CSV file.

        Args:
            source_config: Configuration for the data source, must contain 'benchmark_id'
            data_dir: Directory to download data to

        Returns:
            Path to the downloaded zip file.

        Raises:
            DataSourceError: If download fails
        """
        # Validate config first
        self.validate_config(source_config)

        # Define zip file path
        zip_path = data_dir / "cv_folds_singles_substitutions.zip"
        metadata_path = data_dir / "DMS_substitutions.csv"

        # Check class-level cache first
        if self._cached_zip_path and self._cached_zip_path.exists():
            logger.info(f"Using cached ProteinGym data: {self._cached_zip_path}")
            return self._cached_zip_path

        # Check if file already exists on disk
        if zip_path.exists():
            logger.info(f"ProteinGym data already exists at {zip_path}")
            self._cached_zip_path = zip_path
            return zip_path

        try:
            # Ensure data directory exists
            data_dir.mkdir(parents=True, exist_ok=True)

            logger.info(f"Downloading ProteinGym DMS data from {self.DATASET_URL}...")

            # Download the zip file
            response = requests.get(self.DATASET_URL, stream=True)
            response.raise_for_status()  # Raise exception for bad status codes

            # Save the zip file
            with open(zip_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)

            logger.info(f"Successfully downloaded ProteinGym DMS data to {zip_path}")

            # Cache the path for future use
            self._cached_zip_path = zip_path

            # Download the metadata CSV
            response = requests.get(self.METADATA_URL)
            response.raise_for_status()

            with open(metadata_path, "w", encoding="utf-8") as f:
                f.write(response.text)

            self._cached_metadata = metadata_path

            # Return the path to the zip file (framework will handle extraction)
            return zip_path

        except requests.RequestException as e:
            raise DataSourceError(
                f"Failed to download ProteinGym DMS data: {e}", source_type="proteingym"
            ) from e
        except Exception as e:
            raise DataSourceError(
                f"Unexpected error downloading ProteinGym DMS data: {e}", source_type="proteingym"
            ) from e

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
        if self._cached_leaderboard is not None:
            logger.info("Using cached ProteinGym leaderboard")
            leaderboard_df = self._cached_leaderboard
        else:
            try:
                # Download and parse leaderboard CSV
                leaderboard_df = self._download_and_parse_leaderboard_csv()
                if leaderboard_df is not None and not leaderboard_df.empty:
                    self._cached_leaderboard = leaderboard_df
                else:
                    # If CSV parsing fails, return empty DataFrame
                    logger.warning("CSV parsing failed, returning empty leaderboard")
                    empty_leaderboard = pd.DataFrame(
                        columns=["teamName", "score", "submissionDate"]
                    )
                    self._cached_leaderboard = empty_leaderboard
                    return empty_leaderboard

            except Exception as e:
                logger.error(f"Error downloading/parsing ProteinGym leaderboard: {e}")
                empty_leaderboard = pd.DataFrame(columns=["teamName", "score", "submissionDate"])
                self._cached_leaderboard = empty_leaderboard
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
