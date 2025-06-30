"""OpenProblems batch_integration grading using scIB metrics."""
import anndata as ad
import scib
import numpy as np
import pandas as pd
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

def grade(submission_path: Path, answers_path: Path) -> float:
    """
    Grade batch integration using scIB metrics.
    
    Args:
        submission_path: Path to submitted h5ad file (agent's integration result)
        answers_path: Path to solution h5ad file (ground truth integration)
        
    Returns:
        float: Combined score from multiple scIB metrics (higher is better)
    """
    
    logger.info(f"Loading submission: {submission_path}")
    logger.info(f"Loading solution: {answers_path}")
    
    # Load the AnnData objects
    adata_integrated = ad.read_h5ad(submission_path)
    adata_solution = ad.read_h5ad(answers_path)
    
    logger.info(f"Submission shape: {adata_integrated.shape}")
    logger.info(f"Solution shape: {adata_solution.shape}")
    
    # Ensure both datasets have required metadata
    if 'batch' not in adata_integrated.obs.columns:
        raise ValueError("Submission must contain 'batch' column in obs")
    if 'cell_type' not in adata_integrated.obs.columns:
        raise ValueError("Submission must contain 'cell_type' column in obs")

    # Bio conservation metrics (higher is better)
    logger.info("Computing bio conservation metrics...")
    
    # Calculate scIB metrics
    metrics = {}

    if adata_integrated.X is None:
        adata_integrated.X = adata_integrated.layers['normalized']
    scib.pp.reduce_data(adata_integrated, n_top_genes=2000, batch_key='batch', pca=True)
    adata_integrated = adata_integrated[:, adata_integrated.var['highly_variable']]
    metrics['asw_label'] = scib.me.silhouette(
        adata_integrated,
        group_key='cell_type',
        embed='X_pca'
    )
    metrics['graph_conn'] = scib.me.graph_connectivity(
        adata_integrated, 
        label_key='cell_type'
    )
    logger.info("Computing batch correction metrics...")
    asw_batch_raw = scib.me.silhouette_batch(
        adata_integrated,
        batch_key='batch', 
        group_key='cell_type',
        embed='X_pca'
    )

    # Convert to 0-1 scale (higher = better batch mixing)
    metrics['asw_batch'] = 1 - abs(asw_batch_raw)
        

    metrics['kbet'] = scib.me.kBET(
        adata_integrated,
        batch_key="batch",
        label_key="cell_type",
        type_="embed",
        embed="X_pca",
        scaled=True,
        verbose=True,
    )
    
    # Log individual metrics
    logger.info("Individual metrics:")
    for metric_name, value in metrics.items():
        logger.info(f"  {metric_name}: {value:.4f}")
    
    # Calculate overall score (average of all metrics)
    # All metrics are scaled to 0-1 where higher is better
    valid_metrics = [v for v in metrics.values() if not np.isnan(v) and v is not None]
    
    if not valid_metrics:
        logger.error("No valid metrics computed!")
        return 0.0
    
    overall_score = np.mean(valid_metrics)
    logger.info(f"Overall score: {overall_score:.4f}")
    
    return float(overall_score)


