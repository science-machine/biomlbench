import subprocess
import shutil
from pathlib import Path
import pandas as pd
import anndata as ad
import numpy as np
from biomlbench.utils import get_logger

logger = get_logger(__name__)

# OpenProblems S3 bucket URL for BMMC CITE-seq dataset
OPENPROBLEMS_S3_BASE = "s3://openproblems-data/resources/task_predict_modality/datasets/openproblems_neurips2021/bmmc_cite/normal/log_cp10k/"

def download_from_s3(s3_url: str, local_path: Path) -> None:
    """Download a file from S3 using AWS CLI."""
    try:
        cmd = ["aws", "s3", "cp", s3_url, str(local_path), "--no-sign-request"]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"Successfully downloaded {s3_url} to {local_path}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to download {s3_url}: {e.stderr}")
        raise RuntimeError(f"S3 download failed: {e.stderr}")

def prepare(raw: Path, public: Path, private: Path) -> None:
    """Prepare the predict_modality data from OpenProblems S3 bucket."""
    public.mkdir(parents=True, exist_ok=True)
    private.mkdir(parents=True, exist_ok=True)
    
    # Files to download
    files_to_download = {
        "train_mod1.h5ad": "Training RNA expression data",
        "train_mod2.h5ad": "Training protein abundance data",
        "test_mod1.h5ad": "Test RNA expression data",
        "test_mod2.h5ad": "Test protein abundance data (ground truth)",
    }
    
    # Download all files
    logger.info("Downloading BMMC CITE-seq data from OpenProblems S3 bucket...")
    for filename, description in files_to_download.items():
        logger.info(f"Downloading {filename}: {description}")
        s3_url = f"{OPENPROBLEMS_S3_BASE}{filename}"
        local_path = raw / filename
        download_from_s3(s3_url, local_path)
    
    # Load the data
    logger.info("Loading AnnData objects...")
    train_mod1 = ad.read_h5ad(raw / "train_mod1.h5ad")
    train_mod2 = ad.read_h5ad(raw / "train_mod2.h5ad")
    test_mod1 = ad.read_h5ad(raw / "test_mod1.h5ad")
    test_mod2 = ad.read_h5ad(raw / "test_mod2.h5ad")
    
    # Log dataset information
    logger.info(f"Training data shapes:")
    logger.info(f"  - RNA (mod1): {train_mod1.shape} (cells x genes)")
    logger.info(f"  - Protein (mod2): {train_mod2.shape} (cells x proteins)")
    logger.info(f"Test data shapes:")
    logger.info(f"  - RNA (mod1): {test_mod1.shape} (cells x genes)")
    logger.info(f"  - Protein (mod2): {test_mod2.shape} (cells x proteins)")
    
    # Log data location info
    logger.info("Data is stored in layers['normalized'] (not .X)")
    
    # Verify that training cells match between modalities
    if not all(train_mod1.obs_names == train_mod2.obs_names):
        raise ValueError("Training cell IDs don't match between modalities!")
    
    # Verify that test cells match between modalities
    if not all(test_mod1.obs_names == test_mod2.obs_names):
        raise ValueError("Test cell IDs don't match between modalities!")
    
    # Extract protein names from training data
    protein_names = train_mod2.var_names.tolist()
    logger.info(f"Number of proteins to predict: {len(protein_names)}")
    
    # Create answers.csv with ground truth protein values
    logger.info("Creating answers.csv with ground truth protein values...")
    
    # Get the expression matrix from test_mod2 - it's in the 'normalized' layer
    if 'normalized' in test_mod2.layers:
        test_protein_values = test_mod2.layers['normalized']
        if hasattr(test_protein_values, 'toarray'):
            # Convert sparse matrix to dense
            test_protein_values = test_protein_values.toarray()
    else:
        raise ValueError("test_mod2 doesn't have 'normalized' layer")
    
    # Create DataFrame with cell IDs and protein values
    answers_df = pd.DataFrame(
        test_protein_values,
        index=test_mod2.obs_names,
        columns=protein_names
    )
    answers_df.reset_index(inplace=True)
    answers_df.rename(columns={'index': 'cell_id'}, inplace=True)
    
    # Save answers
    answers_df.to_csv(private / "answers.csv", index=False)
    logger.info(f"Saved answers.csv with shape: {answers_df.shape}")
    
    # Create sample submission with random values
    logger.info("Creating sample submission...")
    
    # Generate random predictions in a reasonable range
    # Look at the range of values in training protein data
    train_protein_values = train_mod2.layers['normalized']
    if hasattr(train_protein_values, 'toarray'):
        train_protein_values = train_protein_values.toarray()
    
    # Calculate statistics for realistic random values
    protein_means = np.mean(train_protein_values, axis=0)
    protein_stds = np.std(train_protein_values, axis=0)
    
    # Generate random predictions
    n_test_cells = test_mod1.shape[0]
    random_predictions = np.random.normal(
        loc=protein_means,
        scale=protein_stds,
        size=(n_test_cells, len(protein_names))
    )
    
    # Create sample submission DataFrame
    sample_submission = pd.DataFrame(
        random_predictions,
        index=test_mod1.obs_names,
        columns=protein_names
    )
    sample_submission.reset_index(inplace=True)
    sample_submission.rename(columns={'index': 'cell_id'}, inplace=True)
    
    # Save sample submission
    sample_submission.to_csv(public / "sample_submission.csv", index=False)
    logger.info(f"Saved sample_submission.csv with shape: {sample_submission.shape}")
    
    # Copy training data and test RNA data to public folder
    logger.info("Copying public data files...")
    shutil.copy(raw / "train_mod1.h5ad", public / "train_mod1.h5ad")
    shutil.copy(raw / "train_mod2.h5ad", public / "train_mod2.h5ad")
    shutil.copy(raw / "test_mod1.h5ad", public / "test_mod1.h5ad")
    
    # Create a metadata file with dataset information
    logger.info("Creating dataset metadata...")
    metadata = {
        "n_train_cells": train_mod1.shape[0],
        "n_test_cells": test_mod1.shape[0],
        "n_genes": train_mod1.shape[1],
        "n_proteins": train_mod2.shape[1],
        "gene_names": train_mod1.var_names.tolist()[:10] + ["..."],  # Sample of gene names
        "protein_names": protein_names,
    }
    
    # Save metadata as a text file for reference
    with open(public / "dataset_info.txt", "w") as f:
        f.write("BMMC CITE-seq Dataset Information\n")
        f.write("=================================\n\n")
        f.write(f"Training cells: {metadata['n_train_cells']:,}\n")
        f.write(f"Test cells: {metadata['n_test_cells']:,}\n")
        f.write(f"Number of genes (RNA): {metadata['n_genes']:,}\n")
        f.write(f"Number of proteins: {metadata['n_proteins']:,}\n")
        f.write(f"\nProtein targets to predict ({metadata['n_proteins']} total):\n")
        for i, protein in enumerate(protein_names):
            f.write(f"  {i+1}. {protein}\n")
        f.write(f"\nData format:\n")
        f.write("- All data is in H5AD (AnnData) format\n")
        f.write("- Expression values are log-transformed and normalized to CP10K\n")
        f.write("- Data is stored in layers['normalized'] (not in .X)\n")
        f.write("- Training data has matched cells between RNA and protein modalities\n")
        f.write("- Task: Predict protein abundance from RNA expression\n")
    
    logger.info("Data preparation complete!")
    logger.info(f"Public files in: {public}")
    logger.info(f"Private files in: {private}")

if __name__ == "__main__":
    # For testing
    prepare(Path("raw"), Path("public"), Path("private")) 