# TODO: This needs to be thoroughly tested
# TODO: We will also need to accept that the agent could peak at half the true results (maybe less comparaible to the public leaderboard)
# TODO: on the other hand, we could basically pretend that the public methods had the same opportunity to peak but they just chose not to because none of them used supervised ML

"""OpenProblems batch_integration data preparation."""
from pathlib import Path
import boto3
import logging
import scanpy as sc
import requests
import pandas as pd
import numpy as np
import collections
import scipy
import anndata
import os
from typing import Union
import pathlib
from tqdm import tqdm
import collections

logger = logging.getLogger(__name__)

URL = "https://figshare.com/ndownloader/files/37593188"

np.random.seed(42)


def map_gene_symbols(adata, map_filename: Union[str, pathlib.Path]):
    """Maps gene symbols from aliases to standard symbols

    Genes that map many-to-one are summed.
    Genes that map one-to-many are duplicated.

    Parameters
    ----------
    adata : anndata.AnnData
    map_filename : PathLike
        Path to csv containing gene symbol map with two columns, `gene` and `alias`

    Returns
    -------
    adata : anndata.AnnData
    """
    map_df = pd.read_csv(map_filename)
    var = adata.var.rename_axis("alias", axis=0)[[]]
    gene_match_idx = np.isin(var.index, map_df["gene"])
    var_gene_match, var = var.loc[gene_match_idx].copy(), var.loc[~gene_match_idx]
    alias_match_idx = np.isin(var.index, map_df["alias"])
    var_alias_match, var_no_map = (
        var.loc[alias_match_idx].copy(),
        var.loc[~alias_match_idx].copy(),
    )

    # fill 'gene' column
    var_alias_match = var_alias_match.reset_index().merge(map_df, on="alias", how="left")
    var_gene_match["gene"] = var_gene_match.index
    var_no_map["gene"] = var_no_map.index

    var_dealiased = pd.concat(
        [var_gene_match.reset_index(), var_no_map.reset_index(), var_alias_match]
    )
    duplicate_idx = var_dealiased["gene"].duplicated(keep=False)
    var_dealiased_many_to_one, var_dealiased_one_to_any = (
        var_dealiased.loc[duplicate_idx],
        var_dealiased.loc[~duplicate_idx],
    )

    adata_one_to_any = adata[:, var_dealiased_one_to_any["alias"]]
    adata_one_to_any.var.index = var_dealiased_one_to_any["gene"]

    many_to_one_genes = var_dealiased_many_to_one["gene"].unique()
    many_to_one_X = []
    many_to_one_layers = collections.defaultdict(list)
    for gene in tqdm(var_dealiased_many_to_one["gene"].unique()):
        gene_aliases = var_dealiased_many_to_one.loc[
            var_dealiased_many_to_one["gene"] == gene, "alias"
        ]
        adata_gene = adata[:, gene_aliases]
        many_to_one_X.append(scipy.sparse.coo_matrix(adata_gene.X.sum(axis=1)))
        for layer_name, layer in adata_gene.layers.items():
            many_to_one_layers[layer_name].append(scipy.sparse.coo_matrix(adata_gene.X.sum(axis=1)))

    return anndata.AnnData(
        X=scipy.sparse.hstack([adata_one_to_any.X] + many_to_one_X).tocsr(),
        obs=adata.obs,
        var=pd.DataFrame(index=np.concatenate([adata_one_to_any.var.index, many_to_one_genes])),
        layers={
            layer_name: scipy.sparse.hstack(
                [adata_one_to_any.layers[layer_name]] + many_to_one_layers[layer_name]
            ).tocsr()
            for layer_name in adata.layers
        },
        uns=adata.uns,
        obsm=adata.obsm,
    )


def download_file(url: str, output_path: Path) -> None:
    """Download a file from a URL to a local path."""
    # Ensure the output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
    # Download the file
    os.system(f"wget {url} -O {output_path}")


def filter_genes_cells(adata):
    """Remove empty cells and genes."""
    import scanpy as sc

    if "var_names_all" not in adata.uns:
        # fill in original var names before filtering
        adata.uns["var_names_all"] = adata.var.index.to_numpy()
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.filter_cells(adata, min_counts=2)


def subsample_even(adata, n_obs, even_obs):
    """Subsample a dataset evenly across an obs.

    Parameters
    ----------
    adata : AnnData
    n_obs : int
        Total number of cells to retain
    even_obs : str
        `adata.obs[even_obs]` to be subsampled evenly across partitions.

    Returns
    -------
    adata : AnnData
        Subsampled AnnData object
    """
    import scanpy as sc

    values = adata.obs[even_obs].unique()
    adatas = []
    n_obs_per_value = n_obs // len(values)
    for v in values:
        adata_subset = adata[adata.obs[even_obs] == v].copy()
        sc.pp.subsample(adata_subset, n_obs=min(n_obs_per_value, adata_subset.shape[0]))
        adatas.append(adata_subset)

    adata_out = anndata.concat(adatas, label="_obs_batch")

    adata_out.uns = adata.uns
    adata_out.varm = adata.varm
    adata_out.varp = adata.varp
    return adata_out


def prepare(raw: Path, prepared: Path) -> None:
    """
    Prepare OpenProblems batch_integration dataset.

    Downloads real h5ad files from S3 for batch integration evaluation.
    - dataset.h5ad: Input data with batch effects (public)
    - solution.h5ad: Ground truth integrated data (private)

    Args:
        raw: Directory for raw data (unused)
        prepared: Directory for prepared data - will create 0/public and 0/private subdirectories
    """
    dataset_dir = prepared / "0"
    public_dir = dataset_dir / "public"
    private_dir = dataset_dir / "private"

    logger.info("Downloading OpenProblems batch_integration h5ad files from S3...")

    logger.info(f"Using dataset: {URL}")

    # Download dataset.h5ad (input with batch effects) from datasets location -> public
    dataset_path = public_dir / "dataset.h5ad"

    logger.info(f"Downloading {URL} to {dataset_path}")
    download_file(URL, dataset_path)

    # Load the dataset
    adata = sc.read_h5ad(dataset_path)

    # Get the map from https://raw.githubusercontent.com/openproblems-bio/openproblems/refs/tags/v1.0.0/openproblems/tasks/_cell_cell_communication/cell_cell_communication_ligand_target/datasets/tnbc_wu2021_gene_symbols.csv
    map_path = public_dir / "tnbc_wu2021_gene_symbols.csv"
    download_file(
        "https://raw.githubusercontent.com/openproblems-bio/openproblems/refs/tags/v1.0.0/openproblems/tasks/_cell_cell_communication/cell_cell_communication_ligand_target/datasets/tnbc_wu2021_gene_symbols.csv",
        map_path,
    )

    adata = map_gene_symbols(adata, map_path)
    adata.uns["merge_keys"] = ["ligand", "target"]

    # Create a public version that is missing the

    # Load the resource from https://raw.githubusercontent.com/saezlab/liana-py/e0c86b15cc2731dde8f4caf63e7918c032f4f2b3/liana/resource/omni_resource.csv
    resource_path = public_dir / "resource.csv"
    download_file(
        "https://raw.githubusercontent.com/saezlab/liana-py/e0c86b15cc2731dde8f4caf63e7918c032f4f2b3/liana/resource/omni_resource.csv",
        resource_path,
    )

    # Load the resource
    resource = pd.read_csv(resource_path, index_col=0)
    resource = resource[resource.resource == "consensus"]
    adata.uns["ligand_receptor_resource"] = resource

    # For the adata.uns['ccc_target'], split it so that 30% of all targets are in the test set
    test_targets = adata.uns["ccc_target"].sample(frac=0.3)
    # For these selected targets, put NaN in the response column
    adata_train = adata.copy()
    adata_train.uns["ccc_target"].loc[test_targets.index, "response"] = np.nan
    adata_test = adata.copy()

    # Save the train and test datasets. Test goes ONLY in private
    adata_train.write_h5ad(public_dir / "dataset.h5ad")
    adata_test.write_h5ad(private_dir / "solution.h5ad")

    # For the adata.uns['ccc_ligand'], split it so that 50% of all ligands are in the test set
    test_ligands = adata.uns["ccc_target"].sample(frac=0.5)

    # Check file size
    dataset_size_mb = dataset_path.stat().st_size / (1024 * 1024)
    logger.info(f"Downloaded dataset.h5ad ({dataset_size_mb:.1f} MB)")

    # Download solution.h5ad (ground truth integration) from task location -> private
    solution_key = f"resources/task_cell_cell_communication/datasets/{dataset_id}/{subset_id}/{normalization_id}/solution.h5ad"
    solution_path = private_dir / "solution.h5ad"

    logger.info(f"Downloading {solution_key} to {solution_path}")
    s3_client.download_file(bucket_name, solution_key, str(solution_path))

    # Check file size
    solution_size_mb = solution_path.stat().st_size / (1024 * 1024)
    logger.info(f"Downloaded solution.h5ad ({solution_size_mb:.1f} MB)")

    # Download and save state.yaml for reference
    try:
        state_key = f"resources/task_cell_cell_communication/datasets/{dataset_id}/{subset_id}/{normalization_id}/state.yaml"
        response = s3_client.get_object(Bucket=bucket_name, Key=state_key)
        state_content = response["Body"].read().decode("utf-8")

        (public_dir / "dataset_info.yaml").write_text(state_content)
        logger.info("Downloaded dataset metadata (state.yaml)")
    except Exception as e:
        logger.warning(f"Could not download state.yaml: {e}")

    # Create task description
    description_content = f"""# OpenProblems Batch Integration

This task evaluates methods for removing batch effects from single-cell data while preserving biological variation.

## Dataset Information

- **Source**: OpenProblems (cell_cell_communication)
- **Input**: `dataset.h5ad` - Single-cell data with ligand-receptor interactions but with test interactions masked (adata.uns['ccc_target'] has NaN in the response column)
- **Expected Output**: `solution.h5ad` - Single-cell data with ligand-receptor interactions with probabilities for each interaction in the adata.uns['ccc_target'].score column
- **Evaluation**: AU-PRC for the scores vs the true interaction labels for the test set. 

## Task Description

Cell-cell communication is a key process in multicellular organisms. This task evaluates methods for predicting ligand-receptor interactions between cells.

## File Sizes

- Dataset: {dataset_size_mb:.1f} MB
- Solution: {solution_size_mb:.1f} MB

Visit https://openproblems.bio for more information about this task.
"""

    (public_dir / "description.md").write_text(description_content)

    logger.info("✅ OpenProblems cell-cell communication data prepared successfully")
    logger.info(f"Dataset: {dataset_path} ({dataset_size_mb:.1f} MB)")
    logger.info(f"Solution: {solution_path} ({solution_size_mb:.1f} MB)")
