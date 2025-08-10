import shutil
import subprocess
from pathlib import Path
from typing import List

import pandas as pd
import scanpy as sc

from biomlbench.utils import get_logger

logger = get_logger(__name__)

# OpenProblems S3 bucket URLs
OPENPROBLEMS_S3_BASE = "s3://openproblems-data/resources/task_perturbation_prediction/datasets/neurips-2023-data/"

def download_from_s3(s3_url: str, local_path: Path) -> None:
    """Download a file from S3 using AWS CLI."""
    try:
        cmd = ["aws", "s3", "cp", s3_url, str(local_path)]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"Successfully downloaded {s3_url} to {local_path}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to download {s3_url}: {e.stderr}")
        raise RuntimeError(f"S3 download failed: {e.stderr}")


def h5ad_to_dataframe(h5ad_path: Path) -> pd.DataFrame:
    """Convert H5AD file to DataFrame with genes as columns."""
    logger.info(f"Loading H5AD file: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    
    # The differential expression data must be in clipped_sign_log10_pval layer
    if not adata.layers or 'clipped_sign_log10_pval' not in adata.layers:
        available_layers = list(adata.layers.keys()) if adata.layers else 'None'
        raise ValueError(f"Required layer 'clipped_sign_log10_pval' not found in {h5ad_path}. Available layers: {available_layers}")
    
    target_layer = 'clipped_sign_log10_pval'
    logger.info(f"Using layer '{target_layer}' for differential expression data")
    
    # Extract the data matrix from the chosen layer
    data_matrix = adata.layers[target_layer]
    
    # Convert sparse matrix to dense if needed
    if hasattr(data_matrix, 'toarray'):
        data_matrix = data_matrix.toarray()
    
    # Create DataFrame with proper structure
    # In AnnData, obs are samples (rows) and var are genes (columns)
    df = pd.DataFrame(
        data_matrix,
        index=adata.obs.index,
        columns=adata.var.index
    )
    
    # Add observation metadata
    for col in adata.obs.columns:
        df[col] = adata.obs[col].values
    
    logger.info(f"Converted H5AD to DataFrame: {df.shape[0]} samples, {df.shape[1]} total columns")
    logger.info(f"Gene columns: {df.shape[1] - len(adata.obs.columns)}, Metadata columns: {len(adata.obs.columns)}")
    return df


def get_gene_columns(df: pd.DataFrame) -> List[str]:
    """Extract gene column names from the dataframe."""
    # Exclude known metadata columns to get just the gene columns
    metadata_cols = {'cell_type', 'sm_name', 'sm_lincs_id', 'SMILES', 'control', 'dose_uM', 'timepoint_hr', 'sm_cell_type', 'split'}
    gene_cols = [col for col in df.columns if col not in metadata_cols]
    return sorted(gene_cols)


def prepare(raw: Path, public: Path, private: Path) -> None:
    """
    Prepare the single-cell perturbations data from OpenProblems S3 bucket.
    
    This downloads the pre-split train/test data from OpenProblems rather than
    creating our own split, since the original Kaggle competition data was incomplete.
    
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
    
    # Download only the essential files we need
    files_to_download = {
        "de_train.h5ad": "Training differential expression data",
        "de_test.h5ad": "Test differential expression data", 
        "id_map.csv": "Test ID mapping",
        "sc_train.h5ad": "Raw single-cell training data",
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
    
    # Load and process the differential expression data
    logger.info("Processing differential expression data...")
    
    # Load training data
    train_df = h5ad_to_dataframe(downloaded_files["de_train.h5ad"])
    logger.info(f"Loaded training data: {len(train_df)} samples")
    
    # Load test data  
    test_df = h5ad_to_dataframe(downloaded_files["de_test.h5ad"])
    logger.info(f"Loaded test data: {len(test_df)} samples")

    
    # Load ID mapping
    id_map = pd.read_csv(downloaded_files["id_map.csv"])
    logger.info(f"Loaded ID mapping: {len(id_map)} test samples")

    
    # Get gene columns
    train_gene_cols = get_gene_columns(train_df)
    test_gene_cols = get_gene_columns(test_df)
    
    # Ensure gene columns match between train and test
    if set(train_gene_cols) != set(test_gene_cols):
        logger.warning("Gene columns differ between train and test data")
        # Use intersection of gene columns
        common_genes = sorted(list(set(train_gene_cols) & set(test_gene_cols)))
        logger.info(f"Using {len(common_genes)} common genes")
        train_gene_cols = common_genes
        test_gene_cols = common_genes
    
    logger.info(f"Working with {len(train_gene_cols)} genes")
    
    # Prepare public files
    logger.info("Creating public files...")
    
    # Save training data as parquet for better performance
    train_df.to_parquet(public / "de_train.parquet", index=False)
    
    # Copy raw single-cell training data for advanced modeling
    if "sc_train.h5ad" in downloaded_files:
        shutil.copy2(downloaded_files["sc_train.h5ad"], public / "sc_train.h5ad")
        logger.info("Copied raw single-cell training data to public directory")
    
    # Save ID mapping
    id_map.to_csv(public / "id_map.csv", index=False)
    
    # Create sample submission with correct format
    sample_submission = pd.DataFrame({
        'id': range(len(test_df)),
        **{gene: 0.0 for gene in test_gene_cols}
    })
    sample_submission.to_csv(public / "sample_submission.csv", index=False)
    
    # Prepare private files
    logger.info("Creating private files...")

    
    # Create answers file from test data
    answers = test_df[test_gene_cols].copy()
    # Use numeric IDs to match sample_submission and id_map
    answers['id'] = range(len(answers))
    
    # Add metadata for verification
    if 'cell_type' in test_df.columns and 'sm_name' in test_df.columns:
        answers['cell_type'] = test_df['cell_type'].values
        answers['sm_name'] = test_df['sm_name'].values

    
    answers.to_csv(private / "answers.csv", index=False)
    
    # Validation checks
    logger.info("Running validation checks...")
    
    # Basic integrity checks
    assert len(train_df) > 0, "Train data is empty"
    assert len(test_df) > 0, "Test data is empty"
    assert len(id_map) == len(test_df), f"ID map length ({len(id_map)}) doesn't match test data ({len(test_df)})"
    assert len(sample_submission) == len(test_df), "Sample submission length doesn't match test data"
    assert len(answers) == len(test_df), "Answers length doesn't match test data"
    
    # Check gene columns
    expected_gene_count = len(test_gene_cols)
    assert len(sample_submission.columns) == expected_gene_count + 1, f"Sample submission should have {expected_gene_count + 1} columns (id + genes)"
    
    # Verify ID mapping matches test data order
    if len(id_map) == len(test_df):
        for i in range(min(10, len(id_map))):  # Check first 10 rows
            if 'cell_type' in test_df.columns and 'sm_name' in test_df.columns:
                assert id_map.iloc[i]['cell_type'] == test_df.iloc[i]['cell_type'], f"Cell type mismatch at row {i}"
                assert id_map.iloc[i]['sm_name'] == test_df.iloc[i]['sm_name'], f"Compound mismatch at row {i}"
    
    # Check data types
    for gene in test_gene_cols[:10]:  # Check first 10 genes
        assert pd.api.types.is_numeric_dtype(answers[gene]), f"Gene {gene} should be numeric in answers"
    
    logger.info("All validation checks passed!")
    
    # Summary statistics
    train_cell_types = train_df['cell_type'].unique() if 'cell_type' in train_df.columns else ['unknown']
    test_cell_types = test_df['cell_type'].unique() if 'cell_type' in test_df.columns else ['unknown']
    
    logger.info(f"Training data: {len(train_df)} samples")
    logger.info(f"Test data: {len(test_df)} samples") 
    logger.info(f"Genes: {len(test_gene_cols)}")
    logger.info(f"Train cell types: {sorted(train_cell_types)}")
    logger.info(f"Test cell types: {sorted(test_cell_types)}")
    
    if 'sm_name' in train_df.columns and 'sm_name' in test_df.columns:
        train_compounds = len(train_df['sm_name'].unique())
        test_compounds = len(test_df['sm_name'].unique())
        logger.info(f"Train compounds: {train_compounds}")
        logger.info(f"Test compounds: {test_compounds}")
    
    logger.info("Data preparation completed successfully using OpenProblems data!") 