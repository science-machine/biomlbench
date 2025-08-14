import subprocess
import shutil
from pathlib import Path
import pandas as pd
import anndata as ad
import numpy as np
from sklearn.model_selection import train_test_split
from biomlbench.utils import get_logger

logger = get_logger(__name__)

# Set a seed
np.random.seed(42)

# OpenProblems S3 bucket URL for Wu et al. breast cancer dataset
OPENPROBLEMS_S3_BASE = "s3://openproblems-data/resources/datasets/openproblems_v1/tnbc_wu2021/log_cp10k/"

def download_from_s3(s3_url: str, local_path: Path) -> None:
    """Download a file from S3 using AWS CLI."""
    try:
        cmd = ["aws", "s3", "cp", s3_url, str(local_path)]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"Successfully downloaded {s3_url} to {local_path}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to download {s3_url}: {e.stderr}")
        raise RuntimeError(f"S3 download failed: {e.stderr}")

def prepare(raw: Path, public: Path, private: Path) -> None:
    """Prepare the cell-cell communication (ligand-target) data from OpenProblems S3 bucket."""
    public.mkdir(parents=True, exist_ok=True)
    private.mkdir(parents=True, exist_ok=True)
    
    # Download TNBC data from S3
    logger.info("Downloading Wu et al. breast cancer data from OpenProblems S3 bucket...")
    dataset_filename = "dataset.h5ad"
    dataset_path = raw / dataset_filename
    
    s3_url = f"{OPENPROBLEMS_S3_BASE}{dataset_filename}"
    download_from_s3(s3_url, dataset_path)
    
    # Load the data
    logger.info("Loading Wu et al. breast cancer dataset...")
    adata = ad.read_h5ad(dataset_path)
    
    # Extract ground truth
    ccc_target = adata.uns["ccc_target"]
    
    # Ensure columns are correct
    if not all(col in ccc_target.columns for col in ["ligand", "target", "response"]):
        raise ValueError("ccc_target must contain columns: ligand, target, response")
    
    # Keep only required columns
    ccc_target = ccc_target[["ligand", "target", "response"]]
    
    # Perform 10/90 train-test split
    logger.info("Performing 10/90 train-test split...")
    # Add index for tracking
    ccc_target = ccc_target.reset_index(drop=True)
    ccc_target["interaction_id"] = range(len(ccc_target))
    
    # Split the data stratified by response to maintain class balance
    train_idx, test_idx = train_test_split(
        ccc_target.index, 
        test_size=0.9, 
        stratify=ccc_target["response"],
        random_state=42
    )
    
    train_data = ccc_target.iloc[train_idx].copy()
    test_data = ccc_target.iloc[test_idx].copy()
    
    logger.info(f"Train set size: {len(train_data)} interactions")
    logger.info(f"Test set size: {len(test_data)} interactions")
    logger.info(f"Train set positive rate: {train_data['response'].mean():.3f}")
    logger.info(f"Test set positive rate: {test_data['response'].mean():.3f}")
    
    # Create answers.csv for grading (only test set)
    logger.info("Creating answers.csv with test set...")
    answers_df = test_data[["ligand", "target", "response", "interaction_id"]].copy()
    answers_df.to_csv(private / "answers.csv", index=False)
    
    # Create sample submission (only for test set)
    logger.info("Creating sample submission for test set...")
    sample_submission = pd.DataFrame({
        "ligand": test_data["ligand"],
        "target": test_data["target"],
        "score": np.random.rand(len(test_data))  # Random scores as example
    })
    sample_submission.to_csv(public / "sample_submission.csv", index=False)
    
    # Prepare public data with training labels
    logger.info("Preparing public data with training labels...")
    adata_public = adata.copy()
    
    # Remove the original full ccc_target
    if "ccc_target" in adata_public.uns:
        del adata_public.uns["ccc_target"]
    if "ccc_pred" in adata_public.uns:
        del adata_public.uns["ccc_pred"]
    
    # Add training data to public adata
    adata_public.uns["ccc_train"] = train_data[["ligand", "target", "response"]]
    
    # Add test ligand-target pairs (without responses) so agents know what to predict
    adata_public.uns["ccc_test_pairs"] = test_data[["ligand", "target"]]
    
    # Save the public data
    adata_public.write_h5ad(public / "tnbc_data.h5ad")
    
    # Also save gene symbol mapping if available
    if "gene_symbols" in adata.var.columns:
        gene_mapping = adata.var[["gene_symbols"]].copy()
        gene_mapping.to_csv(public / "gene_symbols.csv")
    
    # Copy the ligand-receptor prior knowledge resource
    logger.info("Copying ligand-receptor prior knowledge resource...")
    lr_resource_path = Path(__file__).parent / "omnipath.csv.gz"
    if lr_resource_path.exists():
        shutil.copy(lr_resource_path, public / "ligand_receptor_resource.csv.gz")
        logger.info("Ligand-receptor resource copied successfully")
    else:
        logger.warning(f"Ligand-receptor resource not found at {lr_resource_path}")
    
    # Save metadata
    metadata = {
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "cell_types": list(adata.obs["label"].unique()) if "label" in adata.obs else [],
        "n_total_interactions": len(ccc_target),
        "n_train_interactions": len(train_data),
        "n_test_interactions": len(test_data),
        "n_positive_train": int(train_data["response"].sum()),
        "n_negative_train": int((1 - train_data["response"]).sum()),
        "n_positive_test": int(test_data["response"].sum()),
        "n_negative_test": int((1 - test_data["response"]).sum())
    }
    
    logger.info(f"Dataset statistics:")
    logger.info(f"  - Cells: {metadata['n_cells']}")
    logger.info(f"  - Genes: {metadata['n_genes']}")
    logger.info(f"  - Cell types: {len(metadata['cell_types'])}")
    logger.info(f"  - Total interactions: {metadata['n_total_interactions']}")
    logger.info(f"  - Training interactions: {metadata['n_train_interactions']}")
    logger.info(f"  - Test interactions: {metadata['n_test_interactions']}")
    logger.info(f"  - Positive training interactions: {metadata['n_positive_train']}")
    logger.info(f"  - Negative training interactions: {metadata['n_negative_train']}")
    
    # Save metadata as description for agents
    with open(public / "dataset_info.txt", "w") as f:
        f.write("TNBC Cell-Cell Communication Dataset\n")
        f.write("===================================\n\n")
        f.write(f"Number of cells: {metadata['n_cells']:,}\n")
        f.write(f"Number of genes: {metadata['n_genes']:,}\n")
        f.write(f"Number of cell types: {len(metadata['cell_types'])}\n")
        f.write(f"\nData Split:\n")
        f.write(f"- Total interactions: {metadata['n_total_interactions']:,}\n")
        f.write(f"- Training interactions: {metadata['n_train_interactions']:,} (50%)\n")
        f.write(f"- Test interactions: {metadata['n_test_interactions']:,} (50%)\n")
        f.write(f"\nTraining Set:\n")
        f.write(f"- Positive interactions: {metadata['n_positive_train']:,}\n")
        f.write(f"- Negative interactions: {metadata['n_negative_train']:,}\n")
        f.write(f"- Positive rate: {metadata['n_positive_train']/metadata['n_train_interactions']:.3f}\n")
        f.write(f"\nTest Set (to predict):\n")
        f.write(f"- Number of interactions: {metadata['n_test_interactions']:,}\n")
        f.write(f"\nPrior Knowledge:\n")
        f.write(f"- Ligand-receptor resource: ligand_receptor_resource.csv.gz\n")
        f.write(f"- Source: OmniPath consensus database (via LIANA)\n")
        f.write(f"- Contains known ligand-receptor interactions from multiple databases\n")
        f.write(f"\nCell types:\n")
        for ct in sorted(metadata['cell_types']):
            f.write(f"  - {ct}\n")
    
    logger.info("Data preparation complete!")

if __name__ == "__main__":
    # For testing
    prepare(Path("raw"), Path("public"), Path("private")) 