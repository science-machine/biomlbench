"""
Kaggle data source implementation.

This module handles downloading datasets and leaderboards from Kaggle competitions.
"""

import webbrowser
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd
from tenacity import retry, retry_if_exception, stop_after_attempt, wait_fixed

from biomlbench.utils import authenticate_kaggle_api, get_logger
from .base import DataSource, DataSourceError, DataSourceConfigError
from .factory import register_data_source

logger = get_logger(__name__)


def _is_api_exception(exception: Exception) -> bool:
    """Check if exception is a Kaggle API exception."""
    try:
        # only import when necessary; otherwise kaggle asks for API key on import
        from kaggle.rest import ApiException
        return isinstance(exception, ApiException)
    except ImportError:
        return False


def _need_to_accept_rules(error_msg: str) -> bool:
    """Check if error message indicates competition rules need to be accepted."""
    return "You must accept this competition" in error_msg


def _prompt_user_to_accept_rules(competition_id: str) -> None:
    """Prompt user to accept competition rules."""
    response = input("Would you like to open the competition page in your browser now? (y/n): ")

    if response.lower() != "y":
        raise RuntimeError("You must accept the competition rules before downloading the dataset.")

    webbrowser.open(f"https://www.kaggle.com/c/{competition_id}/rules")
    input("Press Enter to continue after you have accepted the rules...")


@register_data_source("kaggle")
class KaggleDataSource(DataSource):
    """
    Data source for Kaggle competitions.
    
    Downloads competition data and leaderboards using the Kaggle API.
    """
    
    def validate_config(self, source_config: Dict[str, Any]) -> bool:
        """
        Validate Kaggle source configuration.
        
        Args:
            source_config: Should contain 'competition_id' key
            
        Returns:
            True if valid
            
        Raises:
            DataSourceConfigError: If configuration is invalid
        """
        if 'competition_id' not in source_config:
            raise DataSourceConfigError(
                "Kaggle data source requires 'competition_id' in configuration",
                source_type="kaggle"
            )
        
        competition_id = source_config['competition_id']
        if not isinstance(competition_id, str) or not competition_id.strip():
            raise DataSourceConfigError(
                "Kaggle 'competition_id' must be a non-empty string",
                source_type="kaggle"
            )
        
        return True
    
    @retry(
        retry=retry_if_exception(_is_api_exception),
        stop=stop_after_attempt(3),  # stop after 3 attempts
        wait=wait_fixed(5),  # wait 5 seconds between attempts
        reraise=True,
    )
    def download(self, source_config: Dict[str, Any], data_dir: Path) -> Optional[Path]:
        """
        Download competition data from Kaggle.
        
        Args:
            source_config: Must contain 'competition_id'
            data_dir: Directory to download data to
            
        Returns:
            Path to downloaded zip file
            
        Raises:
            DataSourceError: If download fails
        """
        self.validate_config(source_config)
        
        competition_id = source_config['competition_id']
        quiet = source_config.get('quiet', False)
        force = source_config.get('force', False)
        
        try:
            if not data_dir.exists():
                data_dir.mkdir(parents=True)

            logger.info(f"Downloading Kaggle dataset for `{competition_id}` to `{data_dir}`...")

            api = authenticate_kaggle_api()

            # only import when necessary; otherwise kaggle asks for API key on import
            from kaggle.rest import ApiException

            try:
                api.competition_download_files(
                    competition=competition_id,
                    path=data_dir,
                    quiet=quiet,
                    force=force,
                )
            except ApiException as e:
                if _need_to_accept_rules(str(e)):
                    logger.warning("You must accept the competition rules before downloading the dataset.")
                    _prompt_user_to_accept_rules(competition_id)
                    # Retry download after accepting rules
                    api.competition_download_files(
                        competition=competition_id,
                        path=data_dir,
                        quiet=quiet,
                        force=force,
                    )
                else:
                    raise e

            zip_files = list(data_dir.glob("*.zip"))

            if len(zip_files) != 1:
                raise DataSourceError(
                    f"Expected to download a single zip file, but found {len(zip_files)} zip files.",
                    source_type="kaggle"
                )

            zip_file = zip_files[0]
            logger.info(f"Successfully downloaded Kaggle dataset: {zip_file}")
            return zip_file

        except Exception as e:
            if isinstance(e, DataSourceError):
                raise
            raise DataSourceError(
                f"Failed to download Kaggle competition '{competition_id}': {e}",
                source_type="kaggle"
            ) from e
    
    def get_leaderboard(self, source_config: Dict[str, Any]) -> pd.DataFrame:
        """
        Get leaderboard from Kaggle competition.
        
        Args:
            source_config: Must contain 'competition_id'
            
        Returns:
            DataFrame with leaderboard data
            
        Raises:
            DataSourceError: If leaderboard cannot be retrieved
        """
        self.validate_config(source_config)
        
        competition_id = source_config['competition_id']
        
        try:
            api = authenticate_kaggle_api()
            leaderboard = api.competition_leaderboard_view(competition=competition_id)
            
            if leaderboard:
                leaderboard = [row.__dict__ for row in leaderboard]
                leaderboard_df = pd.DataFrame(leaderboard)
                # Clean up column names to match expected format
                if "teamNameNullable" in leaderboard_df.columns:
                    leaderboard_df.drop(columns=["teamNameNullable"], inplace=True)
                if "teamName" in leaderboard_df.columns:
                    leaderboard_df.rename(columns={"teamName": "teamName"}, inplace=True)
                
                logger.info(f"Retrieved leaderboard for Kaggle competition '{competition_id}'")
                return leaderboard_df
            else:
                raise DataSourceError(
                    f"Failed to retrieve leaderboard for Kaggle competition '{competition_id}'",
                    source_type="kaggle"
                )
                
        except Exception as e:
            if isinstance(e, DataSourceError):
                raise
            raise DataSourceError(
                f"Failed to get leaderboard for Kaggle competition '{competition_id}': {e}",
                source_type="kaggle"
            ) from e 