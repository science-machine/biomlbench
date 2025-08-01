"""
Data sources package for BioML-bench.

This package provides pluggable data source implementations for different
platforms including Kaggle, Polaris, HTTP, and local files.

The factory pattern allows easy registration and creation of data sources
based on configuration.
"""

from .base import DataSource, DataSourceError, DataSourceConfigError, DataSourceNotFoundError
from .factory import DataSourceFactory, register_data_source

# Import all data source implementations to trigger registration
from .kaggle import KaggleDataSource
from .polaris import PolarisDataSource
from .openproblems import OpenProblemsDataSource

__all__ = [
    # Base classes and interfaces
    "DataSource",
    "DataSourceError", 
    "DataSourceConfigError",
    "DataSourceNotFoundError",
    # Factory
    "DataSourceFactory",
    "register_data_source",
    # Concrete implementations
    "KaggleDataSource",
    "PolarisDataSource",
    "OpenProblemsDataSource",
]


def list_available_sources() -> list[str]:
    """
    Convenience function to list all available data source types.
    
    Returns:
        List of registered data source type strings
    """
    return DataSourceFactory.list_available()


def create_data_source(source_type: str) -> DataSource:
    """
    Convenience function to create a data source instance.
    
    Args:
        source_type: Type of data source to create
        
    Returns:
        Configured data source instance
        
    Raises:
        DataSourceNotFoundError: If source type is not registered
    """
    return DataSourceFactory.create(source_type) 