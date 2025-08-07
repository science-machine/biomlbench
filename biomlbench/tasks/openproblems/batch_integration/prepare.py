"""OpenProblems batch_integration data preparation."""
from pathlib import Path
import boto3
import logging

logger = logging.getLogger(__name__)

def prepare(raw: Path, public: Path, private: Path) -> None:
    """
    Prepare OpenProblems batch_integration dataset.
    
    Downloads real h5ad files from S3 for batch integration evaluation.
    - dataset.h5ad: Input data with batch effects (public)
    - solution.h5ad: Ground truth integrated data (private) 
    
    Args:
        raw: Directory for raw data (unused)
        public: Directory for public data (dataset.h5ad)
        private: Directory for private data (solution.h5ad)
    """
    
    logger.info("Downloading OpenProblems batch_integration h5ad files from S3...")
    
    # Set up S3 client for public bucket (no credentials needed)
    from botocore import UNSIGNED
    from botocore.config import Config
    
    s3_client = boto3.client('s3', 
                            region_name='us-east-1',
                            config=Config(signature_version=UNSIGNED))
    bucket_name = "openproblems-data"
    
    # Use mouse_pancreas_atlas dataset - pick the first normalization method available
    dataset_id = "cellxgene_census"
    subset_id = "mouse_pancreas_atlas" 
    
    # List available normalizations from the datasets location (for input data)
    datasets_prefix = f"resources/datasets/{dataset_id}/{subset_id}/"
    response = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=datasets_prefix, Delimiter='/')
    
    normalizations = []
    for prefix_info in response.get('CommonPrefixes', []):
        norm_path = prefix_info['Prefix']
        norm_id = norm_path.rstrip('/').split('/')[-1]
        if norm_id:  # Skip empty entries
            normalizations.append(norm_id)
    
    if not normalizations:
        raise ValueError(f"No normalizations found for {dataset_id}/{subset_id}")
    
    # Use the first available normalization
    normalization_id = normalizations[0]
    
    logger.info(f"Using dataset: {dataset_id}/{subset_id}/{normalization_id}")
    
    # Download dataset.h5ad (input with batch effects) from datasets location -> public
    dataset_key = f"resources/datasets/{dataset_id}/{subset_id}/{normalization_id}/dataset.h5ad"
    dataset_path = public / "dataset.h5ad"
    
    logger.info(f"Downloading {dataset_key} to {dataset_path}")
    s3_client.download_file(bucket_name, dataset_key, str(dataset_path))
    
    # Check file size
    dataset_size_mb = dataset_path.stat().st_size / (1024 * 1024)
    logger.info(f"Downloaded dataset.h5ad ({dataset_size_mb:.1f} MB)")
    
    # Download solution.h5ad (ground truth integration) from task location -> private
    solution_key = f"resources/task_batch_integration/datasets/{dataset_id}/{subset_id}/{normalization_id}/solution.h5ad"
    solution_path = private / "solution.h5ad"
    
    logger.info(f"Downloading {solution_key} to {solution_path}")
    s3_client.download_file(bucket_name, solution_key, str(solution_path))
    
    # Check file size
    solution_size_mb = solution_path.stat().st_size / (1024 * 1024)
    logger.info(f"Downloaded solution.h5ad ({solution_size_mb:.1f} MB)")
    
    # Download and save state.yaml for reference
    try:
        state_key = f"resources/task_batch_integration/datasets/{dataset_id}/{subset_id}/{normalization_id}/state.yaml"
        response = s3_client.get_object(Bucket=bucket_name, Key=state_key)
        state_content = response['Body'].read().decode('utf-8')
        
        (public / 'dataset_info.yaml').write_text(state_content)
        logger.info("Downloaded dataset metadata (state.yaml)")
    except Exception as e:
        logger.warning(f"Could not download state.yaml: {e}")
    
    # Create task description 
    description_content = f"""# OpenProblems Batch Integration

This task evaluates methods for removing batch effects from single-cell data while preserving biological variation.

## Dataset Information

- **Source**: OpenProblems ({dataset_id}/{subset_id}/{normalization_id})
- **Input**: `dataset.h5ad` - Single-cell data with batch effects
- **Expected Output**: Integrated h5ad file with batch effects removed
- **Evaluation**: Compared against `solution.h5ad` using scIB metrics

## Task Description

Batch integration aims to remove technical differences between batches while preserving biological differences between cell types. Methods are evaluated using multiple metrics including:

- KBET (batch mixing)
- ASW (silhouette width) 
- Graph connectivity
- And other scIB metrics

## File Sizes

- Dataset: {dataset_size_mb:.1f} MB
- Solution: {solution_size_mb:.1f} MB

Visit https://openproblems.bio for more information about this task.
"""
    
    (public / 'description.md').write_text(description_content)
    
    logger.info("✅ OpenProblems batch integration data prepared successfully")
    logger.info(f"Dataset: {dataset_path} ({dataset_size_mb:.1f} MB)")
    logger.info(f"Solution: {solution_path} ({solution_size_mb:.1f} MB)")
