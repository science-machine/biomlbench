# Data Sources API

The data sources module provides a pluggable system for downloading datasets from various biomedical data repositories.

## Base Classes

::: biomlbench.data_sources.base.DataSource

::: biomlbench.data_sources.base.DataSourceError

## Available Data Sources

### Kaggle Data Source

::: biomlbench.data_sources.kaggle.KaggleDataSource

### Polaris Data Source

::: biomlbench.data_sources.polaris.PolarisDataSource

## Factory Functions

::: biomlbench.data_sources.create_data_source

::: biomlbench.data_sources.list_available_sources

## Usage Examples

```python
from biomlbench.data_sources import create_data_source

# Get Kaggle data source
kaggle_source = create_data_source("kaggle")

# Download competition data
kaggle_source.download(
    source_config={"competition_id": "histopathologic-cancer-detection"},
    data_dir=Path("/data")
)

# Get leaderboard
leaderboard = kaggle_source.get_leaderboard(
    source_config={"competition_id": "histopathologic-cancer-detection"}
)
``` 