"""
Abstract base class for data sources in BioML-bench.

This module defines the interface that all data sources must implement,
including methods for downloading data, getting leaderboards, and extracting
human baseline performance data.
"""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, Any, Optional
import pandas as pd


class DataSourceError(Exception):
    """Exception raised for data source related errors."""
    
    def __init__(self, message: str, source_type: str = None):
        self.message = message
        self.source_type = source_type
        super().__init__(self.message)


class DataSource(ABC):
    """Abstract base class for data sources."""
    
    @abstractmethod
    def download(self, source_config: Dict[str, Any], data_dir: Path) -> Optional[Path]:
        """
        Download data from the source.
        
        Args:
            source_config: Configuration specific to this data source
            data_dir: Directory to download data to
            
        Returns:
            Path to downloaded data (or None if no single file)
            
        Raises:
            DataSourceError: If download fails
        """
        pass
    
    @abstractmethod
    def get_leaderboard(self, source_config: Dict[str, Any]) -> pd.DataFrame:
        """
        Get the public leaderboard for the benchmark.
        
        Args:
            source_config: Configuration specific to this data source
            
        Returns:
            DataFrame with columns: teamName, score, submissionDate
            
        Raises:
            DataSourceError: If leaderboard cannot be retrieved
        """
        pass
    
    def get_human_baselines(self, source_config: Dict[str, Any]) -> Optional[pd.DataFrame]:
        """
        Extract human baseline performance data.
        
        Args:
            source_config: Configuration specific to this data source
            
        Returns:
            DataFrame with columns: team_name, score, human_type, source
            Returns None if no human baselines available
            
        Raises:
            DataSourceError: If human baseline extraction fails
        """
        # Default implementation - no human baselines available
        return None
    
    def supports_human_baselines(self) -> bool:
        """Check if this data source can provide human baseline data."""
        return False


class DataSourceNotFoundError(DataSourceError):
    """Exception raised when a data source type is not found."""
    pass


class DataSourceConfigError(DataSourceError):
    """Exception raised when data source configuration is invalid."""
    pass 