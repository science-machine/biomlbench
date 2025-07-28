"""
Kaggle data source implementation.

This module handles downloading datasets and leaderboards from Kaggle competitions.
"""

import webbrowser
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd
from tenacity import retry, retry_if_exception, stop_after_attempt, wait_fixed
from requests.exceptions import HTTPError

from biomlbench.utils import authenticate_kaggle_api, get_logger
from biomlbench.data_sources.base import DataSource, DataSourceError, DataSourceConfigError
from .factory import register_data_source

logger = get_logger(__name__)


def _is_http_exception(exception: Exception) -> bool:
    """Check if exception is an HTTP exception."""
    return isinstance(exception, (HTTPError, ConnectionError))


def _need_to_accept_rules(error_msg: str) -> bool:
    """Check if error message indicates competition rules need to be accepted."""
    return "You must accept this competition" in error_msg or "403" in str(error_msg)


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
            source_config: Should contain 'benchmark_id' key
            
        Returns:
            True if valid
            
        Raises:
            DataSourceConfigError: If configuration is invalid
        """
        if 'benchmark_id' not in source_config:
            raise DataSourceConfigError(
                "Kaggle data source requires 'benchmark_id' in configuration",
                source_type="kaggle"
            )
        
        competition_id = source_config['benchmark_id']
        if not isinstance(competition_id, str) or not competition_id.strip():
            raise DataSourceConfigError(
                "Kaggle 'benchmark_id' must be a non-empty string",
                source_type="kaggle"
            )
        
        return True
    
    @retry(
        retry=retry_if_exception(_is_http_exception),
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
        
        competition_id = source_config['benchmark_id']
        quiet = source_config.get('quiet', False)
        force = source_config.get('force', False)
        
        try:
            if not data_dir.exists():
                data_dir.mkdir(parents=True)

            logger.info(f"Downloading Kaggle dataset for `{competition_id}` to `{data_dir}`...")

            api = authenticate_kaggle_api()

            try:
                api.competition_download_files(
                    competition=competition_id,
                    path=data_dir,
                    quiet=quiet,
                    force=force,
                )
            except Exception as e:
                error_msg = str(e)
                if _need_to_accept_rules(error_msg):
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
            source_config: Must contain 'benchmark_id'
            
        Returns:
            DataFrame with leaderboard data
            
        Raises:
            DataSourceError: If leaderboard cannot be retrieved
        """
        self.validate_config(source_config)
        
        competition_id = source_config['benchmark_id']
        
        try:
            api = authenticate_kaggle_api()
            leaderboard = api.competition_leaderboard_view(competition=competition_id)
            
            if leaderboard:
                leaderboard = [row.__dict__ for row in leaderboard]
                leaderboard_df = pd.DataFrame(leaderboard)
                
                # Normalize column names by removing underscores for consistency
                column_mapping = {}
                for col in leaderboard_df.columns:
                    if col.startswith('_'):
                        new_col = col[1:]  # Remove leading underscore
                        column_mapping[col] = new_col
                
                if column_mapping:
                    leaderboard_df = leaderboard_df.rename(columns=column_mapping)
                
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
    
    def get_human_baselines(self, source_config: Dict[str, Any]) -> Optional[pd.DataFrame]:
        """
        Extract human baselines from Kaggle public leaderboard.
        
        Filters the public leaderboard to identify likely human participants
        and categorizes them by performance level.
        """
        try:
            # Get the full public leaderboard
            leaderboard_df = self.get_leaderboard(source_config)
            
            if leaderboard_df.empty:
                return None
            
            # Filter for likely human teams (heuristic-based)
            human_teams = self._filter_human_teams(leaderboard_df)
            
            if human_teams.empty:
                return None
            
            # Categorize humans by performance level
            human_baselines = self._categorize_human_performance(human_teams)
            
            # Add source metadata
            human_baselines['source'] = 'kaggle_public_leaderboard'
            
            logger.info(f"Extracted {len(human_baselines)} human baseline entries from Kaggle leaderboard")
            return human_baselines
            
        except Exception as e:
            logger.error(f"Failed to extract human baselines from Kaggle: {e}")
            raise e
    
    def _filter_human_teams(self, leaderboard_df: pd.DataFrame) -> pd.DataFrame:
        """Filter leaderboard to identify likely human teams."""
        
        # Create copy to avoid modifying original
        df = leaderboard_df.copy()
        
        # Heuristics to identify human teams vs bots/scripts
        human_indicators = pd.Series(True, index=df.index)
        
        # Filter out obvious bot patterns
        bot_patterns = [
            # Common bot/script patterns in team names
            r'\bbot\b', r'\bapi\b', r'\bscript\b', r'\bauto\b',
            r'\btest\b', r'\bdummy\b', r'\bbaseline\b',
            # Single character or very short names (often bots)
            r'^[a-z]$', r'^[0-9]+$',
            # Very long random-looking strings
            r'^[a-z0-9]{20,}$'
        ]
        
        for pattern in bot_patterns:
            human_indicators &= ~df['team_name'].str.contains(pattern, case=False, na=False, regex=True)
        
        # Filter out teams with suspicious submission patterns
        # (This would require submission frequency data, which Kaggle API may not provide)
        
        # Apply filters
        human_df = df[human_indicators].copy()
        
        logger.debug(f"Filtered {len(df)} total teams to {len(human_df)} likely human teams")
        return human_df
    
    def _categorize_human_performance(self, human_df: pd.DataFrame) -> pd.DataFrame:
        """Categorize human performance into expert/intermediate/novice levels."""
        
        if human_df.empty:
            return pd.DataFrame(columns=['team_name', 'score', 'human_type', 'source'])
        
        # Calculate performance percentiles
        scores = human_df['score'].astype(float)
        
        # Define performance tiers based on percentiles
        expert_threshold = scores.quantile(0.1)      # Top 10%
        intermediate_threshold = scores.quantile(0.5) # Median
        
        # Categorize humans
        human_baselines = []
        
        for _, row in human_df.iterrows():
            score = float(row['score'])
            
            # Determine human type based on performance
            if score <= expert_threshold:
                human_type = 'expert'
            elif score <= intermediate_threshold:
                human_type = 'intermediate'
            else:
                human_type = 'novice'
            
            human_baselines.append({
                'team_name': row['team_name'],
                'score': score,
                'human_type': human_type
            })
        
        return pd.DataFrame(human_baselines)
    
    def supports_human_baselines(self) -> bool:
        """Kaggle supports human baseline extraction from public leaderboards."""
        return True 