"""Manual data source for tasks that handle their own data preparation."""

import pandas as pd
from pathlib import Path
from typing import Dict, Any, Optional

from biomlbench.data_sources.base import DataSource
from biomlbench.data_sources.factory import register_data_source


@register_data_source("manual")
class ManualDataSource(DataSource):
    """
    Data source for manual tasks that handle their own data preparation.
    
    This data source is used for tasks where the data preparation is handled
    entirely by the task's prepare.py script, without needing to download
    data from external sources.
    """
    
    def download(self, source_config: Dict[str, Any], data_dir: Path) -> Optional[Path]:
        """
        Manual tasks don't download data - preparation is handled by prepare.py.
        
        Args:
            source_config: Configuration (not used for manual tasks)
            data_dir: Directory where data would be stored (not used)
            
        Returns:
            None since no download occurs
        """
        # Manual tasks handle data preparation themselves
        return None
    
    def get_leaderboard(self, source_config: Dict[str, Any]) -> pd.DataFrame:
        """
        Manual tasks don't have external leaderboards.
        
        Args:
            source_config: Configuration (not used for manual tasks)
            
        Returns:
            Empty DataFrame with expected leaderboard columns
        """
        # Return empty leaderboard for manual tasks
        return pd.DataFrame(columns=['teamName', 'score', 'submissionDate']) 