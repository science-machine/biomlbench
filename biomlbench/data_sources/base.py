"""
Abstract base class for data sources in BioML-bench.

This module defines the interface that all data sources must implement,
enabling pluggable data acquisition from different platforms.
"""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd


class DataSource(ABC):
    """
    Abstract base class for data sources.
    
    Each data source implementation handles downloading data from a specific 
    platform (Kaggle, Polaris, HTTP, local files, etc.) and provides 
    leaderboard information when available.
    """
    
    @abstractmethod
    def download(self, source_config: Dict[str, Any], data_dir: Path) -> Optional[Path]:
        """
        Download raw data from the source.
        
        Args:
            source_config: Configuration specific to this data source
            data_dir: Directory where raw data should be downloaded
            
        Returns:
            Path to the downloaded data (file or directory), or None if no raw data
            
        Raises:
            DataSourceError: If download fails
        """
        pass
    
    @abstractmethod
    def get_leaderboard(self, source_config: Dict[str, Any]) -> pd.DataFrame:
        """
        Get leaderboard data for the task.
        
        Args:
            source_config: Configuration specific to this data source
            
        Returns:
            DataFrame with leaderboard data (columns: teamName, score, submissionDate)
            
        Raises:
            DataSourceError: If leaderboard cannot be retrieved
        """
        pass
    
    @abstractmethod
    def validate_config(self, source_config: Dict[str, Any]) -> bool:
        """
        Validate that the source configuration is valid for this data source.
        
        Args:
            source_config: Configuration to validate
            
        Returns:
            True if configuration is valid
            
        Raises:
            DataSourceError: If configuration is invalid
        """
        pass
    
    def get_source_type(self) -> str:
        """Get the type identifier for this data source."""
        return self.__class__.__name__.lower().replace('datasource', '')


class DataSourceError(Exception):
    """Exception raised by data source operations."""
    
    def __init__(self, message: str, source_type: str = None):
        self.source_type = source_type
        super().__init__(message)


class DataSourceNotFoundError(DataSourceError):
    """Exception raised when a data source type is not found."""
    pass


class DataSourceConfigError(DataSourceError):
    """Exception raised when data source configuration is invalid."""
    pass 