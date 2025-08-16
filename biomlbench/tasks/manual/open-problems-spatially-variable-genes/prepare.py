import shutil
import subprocess
from pathlib import Path
import pandas as pd
import anndata as ad
import numpy as np

from biomlbench.utils import get_logger

logger = get_logger(__name__)

# OpenProblems S3 bucket URLs for spatially variable genes task
OPENPROBLEMS_S3_BASE = "s3://openproblems-data/resources/task_spatially_variable_genes/datasets/zenodo_spatial/slideseqv2/mouse_cortex/"
OPENPROBLEMS_S3_CEREBELLUM = "s3://openproblems-data/resources/task_spatially_variable_genes/datasets/zenodo_spatial/slideseqv2/mouse_cerebellum/"


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
    """
    Prepare the spatially variable genes data from OpenProblems S3 bucket.
    
    Args:
        raw: Path to raw data directory (not used, data comes from S3)
        public: Path to public data directory (train data and sample submission)
        private: Path to private data directory (ground truth answers)
    """
    # Create output directories
    public.mkdir(parents=True, exist_ok=True)
    private.mkdir(parents=True, exist_ok=True)
    
    # Create a temporary directory for downloads
    temp_dir = Path("temp_svg_download")
    temp_dir.mkdir(exist_ok=True)
    
    # Files to download
    files_to_download = {
        "dataset.h5ad": "Training data with gene expression and spatial coordinates",
        "solution.h5ad": "Ground truth labels for spatially variable genes",
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
    logger.info("Loading AnnData objects...")
    dataset = ad.read_h5ad(downloaded_files["dataset.h5ad"])
    solution = ad.read_h5ad(downloaded_files["solution.h5ad"])
    
    logger.info(f"Dataset shape: {dataset.shape}")
    logger.info(f"Solution shape: {solution.shape}")
    
    # Validate data
    if dataset.n_vars != solution.n_vars:
        raise ValueError(f"Number of genes mismatch: dataset has {dataset.n_vars}, solution has {solution.n_vars}")
    
    # Check for spatial coordinates
    if 'spatial' not in dataset.obsm:
        raise ValueError("No spatial coordinates found in dataset.obsm['spatial']")
    
    # Check for true labels
    if 'true_spatial_var_score' not in solution.var:
        raise ValueError("No true_spatial_var_score found in solution.var")
    
    # Get gene names and extract the scores from the encoded names
    gene_names_with_scores = dataset.var_names.tolist()
    logger.info(f"Processing {len(gene_names_with_scores)} genes")
    
    # Create anonymized gene names (GENE1, GENE2, etc.)
    clean_gene_names = [f"GENE{i+1}" for i in range(len(gene_names_with_scores))]
    
    # Create new AnnData object with anonymized gene names
    logger.info("Creating training data with anonymized gene names...")
    
    # Create new var dataframe with clean names
    new_var = pd.DataFrame(index=clean_gene_names)
    new_var['feature_id'] = clean_gene_names
    new_var['feature_name'] = clean_gene_names
    
    # Create new AnnData with clean gene names
    train_adata = ad.AnnData(
        X=dataset.X,
        obs=dataset.obs.copy(),
        var=new_var,
        uns=dataset.uns.copy(),
        obsm=dataset.obsm.copy()
    )
    
    # Transfer layers
    if 'counts' in dataset.layers:
        train_adata.layers['counts'] = dataset.layers['counts']
    if 'normalized' in dataset.layers:
        train_adata.layers['normalized'] = dataset.layers['normalized']
    
    # Save the anonymized training data
    train_path = public / "train.h5ad"
    train_adata.write_h5ad(train_path)
    logger.info("Created train.h5ad with anonymized gene names")
    
    # Create sample submission with anonymized gene names
    logger.info("Creating sample submission...")
    sample_submission = pd.DataFrame({
        'gene_id': clean_gene_names,
        'spatial_score': np.random.uniform(0, 1, size=len(clean_gene_names))  # Random scores as example
    })
    
    sample_submission_path = public / "sample_submission.csv"
    sample_submission.to_csv(sample_submission_path, index=False)
    logger.info(f"Created sample submission with {len(sample_submission)} genes")
    
    # Create answers file
    logger.info("Creating answers file...")
    
    # Extract true labels from solution
    true_labels = solution.var['true_spatial_var_score'].values
    
    # Check the distribution of true scores
    unique_labels = np.unique(true_labels)
    logger.info(f"Unique true spatial variability scores: {unique_labels}")
    
    answers = pd.DataFrame({
        'gene_id': clean_gene_names,
        'true_spatial_var_score': true_labels
    })
    
    answers_path = private / "answers.csv"
    answers.to_csv(answers_path, index=False)
    logger.info(f"Created answers file with {len(answers)} genes")
    
    # Log score distribution
    n_spatial = (true_labels == 1.0).sum()
    n_non_spatial = (true_labels == 0.0).sum()
    logger.info(f"Score distribution: {n_spatial} highly spatial (1.0), {n_non_spatial} non-spatial (0.0)")
    
    # Validate alignment
    if not all(sample_submission['gene_id'] == answers['gene_id']):
        raise ValueError("Gene ID mismatch between sample submission and answers")
    
    # Download and prepare the cerebellum dataset for testing/validation
    logger.info("Downloading cerebellum dataset for method validation...")
    cerebellum_files = {
        "dataset.h5ad": "Cerebellum training data",
        "solution.h5ad": "Cerebellum ground truth labels",
    }
    
    for filename, description in cerebellum_files.items():
        s3_url = OPENPROBLEMS_S3_CEREBELLUM + filename
        local_path = temp_dir / f"cerebellum_{filename}"
        
        try:
            download_from_s3(s3_url, local_path)
            logger.info(f"Downloaded cerebellum {description}")
            
            # Process files to remove score suffixes from gene names
            if filename == "dataset.h5ad":
                # Load and clean the cerebellum training data
                cerebellum_data = ad.read_h5ad(local_path)
                
                # Create anonymized gene names (GENE1, GENE2, etc.)
                cerebellum_genes_with_scores = cerebellum_data.var_names.tolist()
                cerebellum_clean_genes = [f"GENE{i+1}" for i in range(len(cerebellum_genes_with_scores))]
                
                # Create new var dataframe with clean names
                cerebellum_var = pd.DataFrame(index=cerebellum_clean_genes)
                cerebellum_var['feature_id'] = cerebellum_clean_genes
                cerebellum_var['feature_name'] = cerebellum_clean_genes
                
                # Create new AnnData with clean gene names
                cerebellum_train = ad.AnnData(
                    X=cerebellum_data.X,
                    obs=cerebellum_data.obs.copy(),
                    var=cerebellum_var,
                    uns=cerebellum_data.uns.copy(),
                    obsm=cerebellum_data.obsm.copy()
                )
                
                # Transfer layers
                if 'counts' in cerebellum_data.layers:
                    cerebellum_train.layers['counts'] = cerebellum_data.layers['counts']
                if 'normalized' in cerebellum_data.layers:
                    cerebellum_train.layers['normalized'] = cerebellum_data.layers['normalized']
                
                cerebellum_train.write_h5ad(public / "cerebellum_train.h5ad")
                
            elif filename == "solution.h5ad":
                # For the labels file, we need to keep the scores but clean the gene names
                cerebellum_solution = ad.read_h5ad(local_path)
                
                # Create anonymized gene names (GENE1, GENE2, etc.)
                cerebellum_genes_with_scores = cerebellum_solution.var_names.tolist()
                cerebellum_clean_genes = [f"GENE{i+1}" for i in range(len(cerebellum_genes_with_scores))]
                
                # Create new var dataframe with clean names but keep scores
                cerebellum_solution_var = pd.DataFrame(index=cerebellum_clean_genes)
                cerebellum_solution_var['feature_id'] = cerebellum_clean_genes
                cerebellum_solution_var['feature_name'] = cerebellum_clean_genes
                if 'orig_feature_name' in cerebellum_solution.var.columns:
                    cerebellum_solution_var['orig_feature_name'] = cerebellum_solution.var['orig_feature_name'].values
                cerebellum_solution_var['true_spatial_var_score'] = cerebellum_solution.var['true_spatial_var_score'].values
                
                # Create new AnnData with clean gene names
                cerebellum_labels = ad.AnnData(
                    X=cerebellum_solution.X,
                    obs=cerebellum_solution.obs.copy(),
                    var=cerebellum_solution_var,
                    uns=cerebellum_solution.uns.copy(),
                    obsm=cerebellum_solution.obsm.copy()
                )
                
                cerebellum_labels.write_h5ad(public / "cerebellum_labels.h5ad")
                
        except Exception as e:
            logger.error(f"Failed to download cerebellum {filename}: {e}")
            raise RuntimeError(f"Failed to download cerebellum {filename}: {e}")
    
    logger.info("Created cerebellum validation dataset with labels")
    
    # Clean up temporary directory
    shutil.rmtree(temp_dir)
    logger.info("Cleaned up temporary files")
    
    logger.info("Data preparation complete!")
    
    # Print summary
    print("\nData preparation summary:")
    print(f"- Number of genes: {len(clean_gene_names)}")
    print(f"- Number of spots/cells: {dataset.n_obs}")
    print(f"- Spatial coordinates shape: {dataset.obsm['spatial'].shape}")
    print(f"- Spatially variable genes: {n_spatial} ({n_spatial/len(clean_gene_names)*100:.1f}%)")
    print(f"- Non-spatial genes: {n_non_spatial} ({n_non_spatial/len(clean_gene_names)*100:.1f}%)")
    print("\n- Additional cerebellum dataset provided for method validation with ground truth labels")


if __name__ == "__main__":
    # For testing
    import sys
    if len(sys.argv) != 4:
        print("Usage: python prepare.py <raw_dir> <public_dir> <private_dir>")
        sys.exit(1)
    
    raw = Path(sys.argv[1])
    public = Path(sys.argv[2])
    private = Path(sys.argv[3])
    
    prepare(raw, public, private) 