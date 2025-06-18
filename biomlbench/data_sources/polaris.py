"""
Polaris data source implementation.

This module handles downloading datasets from Polaris Hub and provides
leaderboard information from the Polaris platform.
"""

import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd

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

    def __init__(self):
        self.conda_env = "polarishub"  # Default environment name

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

        # Optional environment override
        if "conda_env" in source_config:
            self.conda_env = source_config["conda_env"]

        return True

    def _run_polaris_script(self, script_content: str) -> str:
        """
        Run a Python script in the polarishub conda environment.

        Args:
            script_content: Python script to execute

        Returns:
            Script output

        Raises:
            DataSourceError: If script execution fails
        """
        with tempfile.NamedTemporaryFile(mode="w", suffix=".py", delete=False) as f:
            f.write(script_content)
            script_path = f.name

        try:
            result = subprocess.run(
                [
                    "bash",
                    "-c",
                    f"source ~/miniconda3/etc/profile.d/conda.sh && "
                    f"conda activate {self.conda_env} && "
                    f"python {script_path}",
                ],
                capture_output=True,
                text=True,
                check=True,
            )

            return result.stdout

        except subprocess.CalledProcessError as e:
            raise DataSourceError(
                f"Polaris script execution failed: {e.stderr}", source_type="polaris"
            ) from e
        finally:
            # Clean up temporary script
            Path(script_path).unlink()

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

        # Create script to download and save Polaris data
        download_script = f"""
import polaris as po
import pandas as pd
from pathlib import Path

try:
    # Load benchmark from Polaris Hub
    benchmark = po.load_benchmark('{benchmark_id}')
    train, _ = benchmark.get_train_test_split()
    # Grab the test set with the targets included
    test = benchmark._get_test_sets(hide_targets=False)['test']

    # Convert to dataframes
    df_train = train.as_dataframe()
    df_test = test.as_dataframe()

    # Save to data directory
    data_dir = Path('{data_dir}')
    df_train.to_csv(data_dir / 'polaris_train_data.csv', index=False)
    df_test.to_csv(data_dir / 'polaris_test_data.csv', index=False)

    print(f"SUCCESS: Downloaded train: {{len(df_train)}} samples")
    print(f"SUCCESS: Downloaded test: {{len(df_test)}} samples")
    print(f"SUCCESS: Saved to {{data_dir}}")
    
except Exception as e:
    print(f"ERROR: {{str(e)}}")
    raise
"""

        try:
            output = self._run_polaris_script(download_script)

            # Check if download was successful
            if "SUCCESS:" not in output:
                raise DataSourceError(f"Polaris download failed: {output}", source_type="polaris")

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
            breakpoint()
            return self._extract_table_data(table, benchmark_id)

        except Exception as e:
            logger.debug(f"Requests approach failed for {benchmark_id}: {e}")
            return None

    def _extract_table_data(self, table, benchmark_id: str) -> Optional[pd.DataFrame]:
        """
        Extract data from a BeautifulSoup table element.

        Args:
            table: BeautifulSoup table element
            benchmark_id: Benchmark ID for context

        Returns:
            DataFrame with leaderboard data or None if extraction fails
        """
        try:
            # Extract headers - only non-empty ones
            header_cells = table.select('thead th[data-slot="table-head"]')
            headers = []
            header_indices = []  # Track which column indices have actual headers

            for i, th in enumerate(header_cells):
                text = th.get_text(strip=True)
                if text and text not in ["", " "]:  # Skip empty headers
                    headers.append(text)
                    header_indices.append(i)

            if len(headers) < 2:  # Need at least rank/name and one metric
                logger.warning(f"Insufficient headers found for {benchmark_id}: {headers}")
                return None

            logger.debug(f"Found headers: {headers} at indices: {header_indices}")

            # Extract data rows
            rows = []
            data_rows = table.select('tbody tr[data-slot="table-row"]')

            for tr in data_rows:
                cells = tr.select('td[data-slot="table-cell"]')

                # Extract only cells that correspond to headers we found
                row_data = []
                for i in header_indices:
                    if i < len(cells):
                        cell = cells[i]
                        # Try different strategies to extract meaningful content
                        text = self._extract_cell_content(cell)
                        row_data.append(text)
                    else:
                        row_data.append("")  # Missing cell

                # Only add rows that have non-empty content in key columns
                if len(row_data) == len(headers) and any(cell.strip() for cell in row_data):
                    row_dict = dict(zip(headers, row_data))
                    rows.append(row_dict)
                    logger.debug(f"Extracted row: {row_dict}")

            if not rows:
                logger.warning(f"No valid data rows found for {benchmark_id}")
                return None

            logger.info(f"Successfully extracted {len(rows)} rows for {benchmark_id}")

            # Convert to standard leaderboard format
            return self._normalize_polaris_leaderboard(pd.DataFrame(rows), benchmark_id)

        except Exception as e:
            logger.warning(f"Failed to extract table data for {benchmark_id}: {e}")
            return None

    def _extract_cell_content(self, cell) -> str:
        """
        Extract meaningful content from a table cell with various strategies.

        Args:
            cell: BeautifulSoup cell element

        Returns:
            Extracted text content
        """
        # Strategy 1: Look for specific content patterns

        # Score/metric values (usually in divs with specific classes)
        score_div = cell.select_one("div.text-center.font-medium")
        if score_div:
            return score_div.get_text(strip=True)

        # Model/team names (usually in <p> tags)
        name_p = cell.select_one("p")
        if name_p:
            return name_p.get_text(strip=True)

        # References (look for links)
        links = cell.select("a[href]")
        if links:
            if len(links) == 1:
                return f"Link: {links[0].get('href', '')}"
            else:
                # Multiple references
                ref_types = []
                for link in links:
                    link_text = link.get_text(strip=True).lower()
                    if "code" in link_text:
                        ref_types.append("Code")
                    elif "paper" in link_text:
                        ref_types.append("Paper")
                    else:
                        ref_types.append("Link")
                return ", ".join(ref_types) if ref_types else "References"

        # Check for "No references" text
        if "no references" in cell.get_text(strip=True).lower():
            return "No references"

        # Strategy 2: Just get all text, but clean it up
        text = cell.get_text(strip=True)

        # Clean up whitespace and empty content
        text = " ".join(text.split())  # Normalize whitespace

        return text if text else ""

    def _normalize_polaris_leaderboard(self, df: pd.DataFrame, benchmark_id: str) -> pd.DataFrame:
        """
        Normalize scraped Polaris leaderboard to standard format.

        Args:
            df: Raw scraped dataframe
            benchmark_id: Benchmark ID for context

        Returns:
            DataFrame with columns: teamName, score, submissionDate
        """
        if df.empty:
            return df

        # Find rank and name columns (usually first few columns)
        rank_col = None
        name_col = None

        for col in df.columns:
            col_lower = col.lower()
            if rank_col is None and ("#" in col or "rank" in col_lower):
                rank_col = col
            elif name_col is None and ("name" in col_lower or "model" in col_lower):
                name_col = col
                break

        if name_col is None:
            # Fallback: assume second column after rank is name
            cols = list(df.columns)
            if len(cols) >= 2:
                name_col = cols[1]

        if name_col is None:
            logger.warning(f"Could not identify name column for {benchmark_id}")
            return pd.DataFrame()  # Return empty

        # Find main metric column (skip rank, name, contributors, references)
        metric_col = None
        skip_patterns = ["#", "rank", "name", "contributor", "reference", "model"]

        for col in df.columns:
            col_lower = col.lower()
            if not any(pattern in col_lower for pattern in skip_patterns):
                # This looks like a metric column
                metric_col = col
                break

        if metric_col is None:
            logger.warning(f"Could not identify metric column for {benchmark_id}")
            # Use placeholder scores
            df["score"] = 0.5
        else:
            # Convert metric to numeric, handling any parsing issues
            df["score"] = pd.to_numeric(df[metric_col], errors="coerce").fillna(0.5)

        # Create normalized leaderboard
        normalized = pd.DataFrame(
            {
                "teamName": df[name_col].astype(str),
                "score": df["score"],
                "submissionDate": "2024-01-01",  # Placeholder date
            }
        )

        # Remove any empty team names
        normalized = normalized[normalized["teamName"].str.strip() != ""]

        logger.info(f"Normalized {len(normalized)} leaderboard entries for {benchmark_id}")
        return normalized
