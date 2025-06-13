"""
Factory for creating data source instances.

This module provides a centralized way to create data source instances
based on configuration, with automatic registration of available sources.
"""

from typing import Dict, Type

from .base import DataSource, DataSourceNotFoundError


class DataSourceFactory:
    """
    Factory for creating data source instances.
    
    Automatically discovers and registers data source classes,
    then creates instances based on source type configuration.
    """
    
    _sources: Dict[str, Type[DataSource]] = {}
    
    @classmethod
    def register(cls, source_type: str, source_class: Type[DataSource]) -> None:
        """
        Register a data source class.
        
        Args:
            source_type: String identifier for the source type
            source_class: Data source class to register
        """
        cls._sources[source_type] = source_class
    
    @classmethod
    def create(cls, source_type: str) -> DataSource:
        """
        Create a data source instance.
        
        Args:
            source_type: Type of data source to create
            
        Returns:
            Configured data source instance
            
        Raises:
            DataSourceNotFoundError: If source type is not registered
        """
        if source_type not in cls._sources:
            available = list(cls._sources.keys())
            raise DataSourceNotFoundError(
                f"Unknown data source type: '{source_type}'. "
                f"Available types: {available}",
                source_type=source_type
            )
        
        source_class = cls._sources[source_type]
        return source_class()
    
    @classmethod
    def list_available(cls) -> list[str]:
        """
        List all available data source types.
        
        Returns:
            List of registered source type strings
        """
        return list(cls._sources.keys())
    
    @classmethod
    def is_available(cls, source_type: str) -> bool:
        """
        Check if a data source type is available.
        
        Args:
            source_type: Source type to check
            
        Returns:
            True if source type is registered
        """
        return source_type in cls._sources


def register_data_source(source_type: str):
    """
    Decorator to automatically register data source classes.
    
    Args:
        source_type: String identifier for the source type
    """
    def decorator(cls: Type[DataSource]) -> Type[DataSource]:
        DataSourceFactory.register(source_type, cls)
        return cls
    return decorator 