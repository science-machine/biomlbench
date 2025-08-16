import shutil
import subprocess
from pathlib import Path
from typing import List, Set

import pandas as pd
import scanpy as sc
import anndata as ad

from biomlbench.utils import get_logger

logger = get_logger(__name__)

# OpenProblems S3 bucket URLs for label projection task
OPENPROBLEMS_S3_BASE = "s3://openproblems-data/resources/task_label_projection/datasets/cellxgene_census/dkd/log_cp10k/"

def download_from_s3(s3_url: str, local_path: Path) -> None:
    """Download a file from S3 using AWS CLI."""
    try:
        cmd = ["aws", "s3", "cp", s3_url, str(local_path), "--no-sign-request"]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"Successfully downloaded {s3_url} to {local_path}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to download {s3_url}: {e.stderr}")
        raise RuntimeError(f"S3 download failed: {e.stderr}")


def get_unique_labels(adata: ad.AnnData) -> List[str]:
    """Extract unique cell type labels from AnnData object."""
    if 'label' not in adata.obs.columns:
        raise ValueError("No 'label' column found in obs")
    
    labels = sorted(adata.obs['label'].unique().tolist())
    return labels


def prepare_test_features(test_adata: ad.AnnData) -> ad.AnnData:
    """Create a version of test data without labels for public use."""
    # Create a copy to avoid modifying original
    public_test = test_adata.copy()
    
    # Remove label column if it exists (it shouldn't in test data)
    if 'label' in public_test.obs.columns:
        logger.warning("Found 'label' column in test data - removing for public version")
        public_test.obs = public_test.obs.drop(columns=['label'])
    
    return public_test


def prepare(raw: Path, public: Path, private: Path) -> None:
    """
    Prepare the label projection data from OpenProblems S3 bucket.
    
    Downloads the DKD (Diabetic Kidney Disease) dataset and prepares it for the
    label projection benchmark task.
    
    Args:
        raw: Path to downloaded raw data (will be used as temp directory)
        public: Path where public data (for agents) should be placed  
        private: Path where private data (answers) should be placed
    """
    public.mkdir(parents=True, exist_ok=True)
    private.mkdir(parents=True, exist_ok=True)
    
    # Use raw directory as temporary download location
    temp_dir = raw
    temp_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("Downloading data from OpenProblems S3 bucket...")
    
    # Files to download for label projection task
    files_to_download = {
        "train.h5ad": "Training data with cell type labels",
        "test.h5ad": "Test data for label prediction",
        "solution.h5ad": "Test data with true labels (private)",
    }
    
    downloaded_files = {}
    
    for filename, description in files_to_download.items():
        s3_url = OPENPROBLEMS_S3_BASE + filename
        local_path = temp_dir / filename
        
        try:
            download_from_s3(s3_url, local_path)
            downloaded_files[filename] = local_path
            logger.info(f"Downloaded {description}")
        except Exception as e:
            logger.error(f"Required file {filename} failed to download: {e}")
            raise RuntimeError(f"Failed to download required file {filename}: {e}")
    
    # Load the data
    logger.info("Loading H5AD files...")
    
    train_adata = ad.read_h5ad(downloaded_files["train.h5ad"])
    logger.info(f"Loaded training data: {train_adata.shape[0]} cells, {train_adata.shape[1]} genes")
    
    test_adata = ad.read_h5ad(downloaded_files["test.h5ad"])
    logger.info(f"Loaded test data: {test_adata.shape[0]} cells, {test_adata.shape[1]} genes")
    
    solution_adata = ad.read_h5ad(downloaded_files["solution.h5ad"])
    logger.info(f"Loaded solution data: {solution_adata.shape[0]} cells")
    
    # Validate data consistency
    logger.info("Validating data consistency...")
    
    # Check that train has labels
    if 'label' not in train_adata.obs.columns:
        raise ValueError("Training data missing 'label' column")
    
    # Check that test doesn't have labels
    if 'label' in test_adata.obs.columns:
        logger.warning("Test data contains labels - this should not happen in production")
    
    # Check that solution has labels
    if 'label' not in solution_adata.obs.columns:
        raise ValueError("Solution data missing 'label' column")
    
    # Verify test and solution have same cells
    if not all(test_adata.obs.index == solution_adata.obs.index):
        raise ValueError("Test and solution data have different cell indices")
    
    # Verify gene sets match
    if not all(train_adata.var.index == test_adata.var.index):
        raise ValueError("Train and test data have different gene sets")
    
    # Get unique labels from training data
    train_labels = get_unique_labels(train_adata)
    logger.info(f"Found {len(train_labels)} unique cell types in training data")
    logger.info(f"Cell types: {train_labels[:10]}..." if len(train_labels) > 10 else f"Cell types: {train_labels}")
    
    # Prepare public files
    logger.info("Creating public files...")
    
    # Save training data as-is (includes labels)
    train_adata.write_h5ad(public / "train.h5ad")
    logger.info("Saved training data to public/train.h5ad")
    
    # Save test data without labels
    public_test = prepare_test_features(test_adata)
    public_test.write_h5ad(public / "test.h5ad")
    logger.info("Saved test features to public/test.h5ad")
    
    # Create sample submission
    sample_submission = pd.DataFrame({
        'cell_id': test_adata.obs.index,
        'label': train_labels[0]  # Use first label as placeholder
    })
    sample_submission.to_csv(public / "sample_submission.csv", index=False)
    logger.info(f"Created sample submission with {len(sample_submission)} cells")
    
    # Create label vocabulary file
    label_vocab = pd.DataFrame({
        'label': train_labels,
        'label_id': range(len(train_labels))
    })
    label_vocab.to_csv(public / "label_vocabulary.csv", index=False)
    logger.info(f"Created label vocabulary with {len(label_vocab)} cell types")
    
    # Prepare private files
    logger.info("Creating private files...")
    
    # Save solution with labels
    solution_adata.write_h5ad(private / "solution.h5ad")
    logger.info("Saved solution data to private/solution.h5ad")
    
    # Create answers CSV for easier evaluation
    answers = pd.DataFrame({
        'cell_id': solution_adata.obs.index,
        'label': solution_adata.obs['label']
    })
    answers.to_csv(private / "answers.csv", index=False)
    logger.info(f"Created answers.csv with {len(answers)} cells")
    
    # Validation checks
    logger.info("Running validation checks...")
    
    # Basic integrity checks
    assert train_adata.shape[0] > 0, "Training data is empty"
    assert test_adata.shape[0] > 0, "Test data is empty"
    assert len(train_labels) > 1, "Need at least 2 cell types for classification"
    assert train_adata.shape[1] == test_adata.shape[1], "Gene counts don't match between train and test"
    
    # Check solution labels are valid
    solution_labels = set(solution_adata.obs['label'].unique())
    train_label_set = set(train_labels)
    unknown_labels = solution_labels - train_label_set
    if unknown_labels:
        logger.warning(f"Solution contains labels not in training set: {unknown_labels}")
    
    # Check data layers
    logger.info("Available data layers:")
    if train_adata.layers:
        logger.info(f"  Training: {list(train_adata.layers.keys())}")
    else:
        logger.info("  Training: No layers, using .X")
    
    if test_adata.layers:
        logger.info(f"  Test: {list(test_adata.layers.keys())}")
    else:
        logger.info("  Test: No layers, using .X")
    
    # Summary statistics
    logger.info("\n=== Data Summary ===")
    logger.info(f"Training cells: {train_adata.shape[0]:,}")
    logger.info(f"Test cells: {test_adata.shape[0]:,}")
    logger.info(f"Number of genes: {train_adata.shape[1]:,}")
    logger.info(f"Number of cell types: {len(train_labels)}")
    
    # Cell type distribution in training
    train_label_counts = train_adata.obs['label'].value_counts()
    logger.info("\nTop 10 cell types in training data:")
    for label, count in train_label_counts.head(10).items():
        logger.info(f"  {label}: {count:,} cells ({count/len(train_adata)*100:.1f}%)")
    
    # Batch information if available
    if 'batch' in train_adata.obs.columns:
        n_batches = train_adata.obs['batch'].nunique()
        logger.info(f"\nNumber of batches in training: {n_batches}")
    
    logger.info("\nData preparation completed successfully!") 