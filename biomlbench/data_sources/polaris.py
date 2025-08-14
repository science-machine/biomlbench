"""
Polaris data source implementation.

This module handles downloading datasets from Polaris Hub and provides
leaderboard information from the Polaris platform.
"""

from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd
import polaris as po

from biomlbench.data_sources.base import DataSource, DataSourceError
from biomlbench.utils import get_logger

from .factory import register_data_source

logger = get_logger(__name__)


@register_data_source("polaris")
class PolarisDataSource(DataSource):
    """
    Data source for Polaris Hub benchmarks.

    Downloads benchmark data and provides leaderboard information
    from the Polaris platform using the polarishub conda environment.
    """

    def validate_config(self, source_config: Dict[str, Any]) -> bool:
        """
        Validate Polaris source configuration.

        Args:
            source_config: Should contain 'benchmark_id' key

        Returns:
            True if valid

        Raises:
            DataSourceConfigError: If configuration is invalid
        """
        if "benchmark_id" not in source_config:
            raise ValueError(
                "Polaris data source requires 'benchmark_id' in configuration",
                source_type="polaris",
            )

        benchmark_id = source_config["benchmark_id"]
        if not isinstance(benchmark_id, str) or not benchmark_id.strip():
            raise ValueError(
                "Polaris 'benchmark_id' must be a non-empty string", source_type="polaris"
            )

        return True

    def download(self, source_config: Dict[str, Any], data_dir: Path) -> Optional[Path]:
        """
        Download benchmark data from Polaris Hub.

        Args:
            source_config: Must contain 'benchmark_id'
            data_dir: Directory to save data to

        Returns:
            Path to the data directory (Polaris doesn't use zip files)

        Raises:
            DataSourceError: If download fails
        """
        self.validate_config(source_config)

        benchmark_id = source_config["benchmark_id"]
        data_dir.mkdir(parents=True, exist_ok=True)

        try:
            benchmark = po.load_benchmark(benchmark_id)
            train, _ = benchmark.get_train_test_split()
            # Grab the test set with the targets included
            test = benchmark._get_test_sets(hide_targets=False)["test"]

            # Convert to dataframes
            df_train = train.as_dataframe()
            df_test = test.as_dataframe()

            # Save to data directory
            df_train.to_csv(data_dir / "polaris_train_data.csv", index=False)
            df_test.to_csv(data_dir / "polaris_test_data.csv", index=False)

            logger.info(f"Successfully downloaded Polaris benchmark '{benchmark_id}' to {data_dir}")

            return data_dir

        except Exception as e:
            if isinstance(e, DataSourceError):
                raise
            raise DataSourceError(
                f"Failed to download Polaris benchmark '{benchmark_id}': {e}", source_type="polaris"
            ) from e

    def get_leaderboard(self, source_config: Dict[str, Any]) -> pd.DataFrame:
        """
        Get leaderboard from Polaris Hub by scraping the website.

        Args:
            source_config: Must contain 'benchmark_id'

        Returns:
            DataFrame with leaderboard data (teamName, score, submissionDate columns)

        Raises:
            DataSourceError: If leaderboard cannot be retrieved
        """
        self.validate_config(source_config)

        benchmark_id = source_config["benchmark_id"]

        try:
            # Try to scrape leaderboard from Polaris Hub website
            leaderboard_df = self._scrape_polaris_leaderboard(benchmark_id)

            if leaderboard_df is not None and not leaderboard_df.empty:
                logger.info(
                    f"Scraped {len(leaderboard_df)} entries from Polaris Hub leaderboard for '{benchmark_id}'"
                )
                return leaderboard_df

            # Fallback to basic structure if scraping fails
            logger.error(f"Could not scrape leaderboard for '{benchmark_id}'")
            raise DataSourceError(
                f"Could not scrape leaderboard for '{benchmark_id}'", source_type="polaris"
            )

        except Exception as e:
            if isinstance(e, DataSourceError):
                raise
            raise DataSourceError(
                f"Failed to get leaderboard for Polaris benchmark '{benchmark_id}': {e}",
                source_type="polaris",
            ) from e

    def _scrape_polaris_leaderboard(self, benchmark_id: str) -> Optional[pd.DataFrame]:
        """
        Scrape leaderboard from Polaris Hub website (first page only).

        Tries lightweight requests + BeautifulSoup first, falls back to Selenium if needed.

        Args:
            benchmark_id: Benchmark ID like 'tdcommons/caco2-wang'

        Returns:
            DataFrame with leaderboard data or None if scraping fails
        """
        # Try lightweight approach first
        result = self._scrape_with_requests(benchmark_id)
        if result is not None:
            logger.info(f"Successfully scraped {benchmark_id} with requests + BeautifulSoup")
            return result

    def _scrape_with_requests(self, benchmark_id: str) -> Optional[pd.DataFrame]:
        """
        Try scraping with lightweight requests + BeautifulSoup.

        Args:
            benchmark_id: Benchmark ID like 'tdcommons/caco2-wang'

        Returns:
            DataFrame with leaderboard data or None if scraping fails
        """
        try:
            import requests
            from bs4 import BeautifulSoup
        except ImportError:
            raise DataSourceError(
                "requests or beautifulsoup4 not available. Install with: pip install requests beautifulsoup4",
                source_type="polaris",
            )

        url = f"https://polarishub.io/benchmarks/{benchmark_id}"

        try:
            # Make HTTP request with proper headers
            headers = {
                "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
            }

            response = requests.get(url, headers=headers, timeout=10)
            response.raise_for_status()

            # Parse HTML
            soup = BeautifulSoup(response.content, "html.parser")

            # Look for the leaderboard table
            table = soup.select_one('table[data-slot="table"]')
            if not table:
                logger.debug(f"No table found with requests for {benchmark_id}")
                return None

            # Check if table has actual data (not just structure)
            data_rows = table.select('tbody tr[data-slot="table-row"]')
            if not data_rows:
                logger.debug(f"Table found but no data rows for {benchmark_id}")
                return None

            # Extract and process the table
            return self._extract_table_data(table, benchmark_id)

        except Exception as e:
            logger.debug(f"Requests approach failed for {benchmark_id}: {e}")
            return None

    def _extract_table_data(self, table, benchmark_id: str) -> Optional[pd.DataFrame]:
        """
        Extract data from a BeautifulSoup table element using pandas read_html.

        Args:
            table: BeautifulSoup table element
            benchmark_id: Benchmark ID for context

        Returns:
            DataFrame with leaderboard data or None if extraction fails
        """
        try:
            # Convert the table element to HTML string and use pandas to parse it
            table_html = str(table)
            # Use pandas read_html to parse the table
            dfs = pd.read_html(table_html, header=0)

            if not dfs:
                logger.warning(f"No tables found in HTML for {benchmark_id}")
                return None

            df = dfs[0]  # Take the first (and likely only) table

            if df.empty:
                logger.warning(f"Empty dataframe extracted for {benchmark_id}")
                return None

            logger.info(f"Successfully extracted table with shape {df.shape} for {benchmark_id}")
            logger.debug(f"Columns: {list(df.columns)}")
            logger.debug(f"First few rows:\n{df.head()}")

            # Convert to standard leaderboard format
            return self._normalize_polaris_leaderboard(df, benchmark_id)

        except Exception as e:
            logger.warning(f"Failed to extract table data for {benchmark_id}: {e}")
            return None

    def _normalize_polaris_leaderboard(self, df: pd.DataFrame, benchmark_id: str) -> pd.DataFrame:
        """
        Normalize scraped Polaris leaderboard to standard format using benchmark's main metric.

        Args:
            df: Raw scraped dataframe
            benchmark_id: Benchmark ID for context

        Returns:
            DataFrame with columns: teamName, score, submissionDate
        """
        if df.empty:
            return df

        # Get the main metric from the benchmark
        try:
            import polaris as po

            benchmark = po.load_benchmark(benchmark_id)
            main_metric_label = benchmark.main_metric.label
            logger.debug(f"Main metric for {benchmark_id}: {main_metric_label}")
        except Exception as e:
            raise DataSourceError(
                f"Could not load benchmark {benchmark_id} to get main metric: {e}",
                source_type="polaris",
            ) from e

        # Find the team/model name column (usually contains "name" or "model")
        name_col = None
        # Remove any columns that match the pattern "Unnamed"
        df = df.loc[:, ~df.columns.str.contains("^Unnamed")]
        for col in df.columns:
            col_lower = str(col).lower()
            if any(keyword in col_lower for keyword in ["name", "model"]):
                name_col = col
                break

        if name_col is None:
            raise DataSourceError(
                f"Could not identify name/model column for {benchmark_id}. "
                f"Available columns: {list(df.columns)}",
                source_type="polaris",
            )

        # Find the score column using the main metric
        score_col = None
        if main_metric_label:
            # Look for exact match first
            for col in df.columns:
                if str(col).strip() == main_metric_label:
                    score_col = col
                    break

            # If no exact match, look for partial match
            if score_col is None:
                for col in df.columns:
                    if main_metric_label.lower() in str(col).lower():
                        score_col = col
                        break

        # If we still don't have a score column, find the first numeric column that's not rank/position
        if score_col is None:
            for col in df.columns:
                col_lower = str(col).lower()
                # Skip obvious non-score columns
                if any(
                    skip in col_lower
                    for skip in [
                        "#",
                        "rank",
                        "position",
                        "name",
                        "model",
                        "reference",
                        "contributor",
                    ]
                ):
                    continue
                # Check if this column contains numeric data
                try:
                    pd.to_numeric(df[col], errors="coerce")
                    score_col = col
                    break
                except:
                    continue

        if score_col is None:
            raise DataSourceError(
                f"Could not identify score column for {benchmark_id}. "
                f"Main metric: '{main_metric_label}', Available columns: {list(df.columns)}",
                source_type="polaris",
            )

        logger.debug(f"Using score column '{score_col}' for {benchmark_id}")
        # Convert to numeric, fail if we can't parse the scores properly
        try:
            scores = pd.to_numeric(df[score_col], errors="raise")
        except (ValueError, TypeError) as e:
            raise DataSourceError(
                f"Could not parse scores from column '{score_col}' for {benchmark_id}: {e}",
                source_type="polaris",
            ) from e

        # Create normalized leaderboard
        team_names = df[name_col].astype(str)

        # Validate that we actually extracted meaningful team names
        valid_names = team_names[
            (team_names.str.strip() != "")
            & (team_names.str.lower() != "nan")
            & (~team_names.isna())
        ]

        if len(valid_names) == 0:
            raise DataSourceError(
                f"No valid team names found for {benchmark_id}. "
                f"All team names are empty, 'nan', or null. "
                f"Raw name column '{name_col}' values: {team_names.tolist()[:10]}",
                source_type="polaris",
            )

        if len(valid_names) < len(team_names) * 0.5:  # More than 50% invalid names
            raise DataSourceError(
                f"Most team names are invalid for {benchmark_id}. "
                f"Only {len(valid_names)}/{len(team_names)} names are valid. "
                f"Raw name column '{name_col}' values: {team_names.tolist()[:10]}",
                source_type="polaris",
            )

        normalized = pd.DataFrame(
            {
                "teamName": team_names,
                "score": scores,
                "submissionDate": "2024-01-01",  # Placeholder date
            }
        )

        # Remove any empty/invalid team names
        normalized = normalized[
            (normalized["teamName"].str.strip() != "")
            & (normalized["teamName"].str.lower() != "nan")
            & (~normalized["teamName"].isna())
        ]

        logger.info(f"Normalized {len(normalized)} leaderboard entries for {benchmark_id}")
        logger.debug(f"Used columns: name='{name_col}', score='{score_col}'")
        return normalized
