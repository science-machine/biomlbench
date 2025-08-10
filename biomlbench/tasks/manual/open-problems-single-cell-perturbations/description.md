# OpenProblems – Single-Cell Perturbation Prediction

## Overview

This task aims to predict how small molecules change gene expression in different cell types. You will predict differential expression values for cell-type and compound combinations not seen during training.

## Task Description

Given training data showing how various compounds affect gene expression in different cell types, predict the differential expression response for new compound-cell type combinations in the test set.

## Dataset

The dataset contains differential expression measurements from a single-cell perturbational experiment in human peripheral blood mononuclear cells (PBMCs).

### Experimental Setup
- **144 compounds** from the LINCS Connectivity Map dataset
- **Treatment duration**: 24 hours
- **Cell types**: B cells and Myeloid cells (test set), plus T cells, NK cells, and regulatory T cells (training set)
- **Three healthy human donors**

### Data Splits
- **Training**: All compounds in T cells, NK cells, regulatory T cells + subset of compounds in B/Myeloid cells
- **Test**: Remaining compounds in B cells and Myeloid cells

## Input Files

You will receive the following files:

### `de_train.parquet`
Training differential expression data with:
- **Rows**: Cell type and compound combinations  
- **Columns**: Gene expression values (5317 genes) + metadata
- **Values**: Differential expression scores (clipped_sign_log10_pval metric)

The training data contains multiple layers of differential expression statistics:
- `clipped_sign_log10_pval`: Primary target (`sign_log10_pval` values clipped between -4 and 4)
- `sign_log10_pval`: Unclipped differential expression values (-log10(p-value) × sign(log_fold_change))
- `logFC`: Log fold change values
- `P.Value`: Raw p-values  
- `adj.P.Value`: Adjusted p-values

### `sc_train.h5ad` (It's not required to use this data but could be helpful for certain approaches.)
Raw single-cell training data in AnnData format for advanced modeling approaches. It was from this data that the differential expression data was extracted (for training samples). adata.X is in raw counts format. adata.obs contains cell type and compount information.  adata.var_names contains the gene names that are also used in the `de_train.parquet` file.

### `id_map.csv`
Maps test set IDs to their corresponding cell types and compounds:

```
id,sm_name,cell_type
0,TIE2 Kinase Inhibitor,B cells
1,TIE2 Kinase Inhibitor,Myeloid cells
...
```

### `sample_submission.csv`
Template showing the required submission format with all gene columns initialized to 0.0.

## Target Variable

Predict the **clipped_sign_log10_pval** values for each gene, calculated as:

```
clipped_sign_log10_pval = clip(−log10(p-value) × sign(log_fold_change), -4, 4)
```

Where:
- Values are clipped to the range [-4, 4] for stability
- Positive values indicate upregulation in treatment vs control
- Negative values indicate downregulation in treatment vs control
- Magnitude reflects statistical significance (higher absolute values = more significant)

## Submission Format

Submit a CSV file with:
- **`id`** column: Test sample IDs (0 to 150)
- **Gene columns**: Predicted differential expression values for 5317 genes
- **Example**: `id,AAK1,AAMP,AASDH,...,ZZEF1`

## Evaluation

Performance is measured using **Mean Rowwise Root Mean Squared Error (MRRMSE)**:
1. Calculate RMSE for each gene across all test samples
2. Take the mean of these gene-wise RMSE values
3. Lower scores are better

Here is the code for the metric:

```python

def mean_rowwise_rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:

    # Calculate squared differences
    squared_diff = (y_true - y_pred) ** 2
    
    # Calculate mean squared error for each row (across genes)
    mse_per_row = np.mean(squared_diff, axis=1)
    
    # Calculate RMSE for each row
    rmse_per_row = np.sqrt(mse_per_row)
    
    # Calculate mean across all rows
    mrrmse = np.mean(rmse_per_row)
    
    return mrrmse

# Prepare data for evaluation
gene_cols = ['AAK1', 'AAMP', 'AASDH', ..., 'ZYX', 'ZZEF1']
y_true = answers[gene_cols].values
y_pred = submission[gene_cols].values

# Calculate MRRMSE
score = mean_rowwise_rmse(y_true, y_pred)

```

This metric treats all genes equally regardless of expression level.

## Scientific Context

This task addresses a key challenge in drug discovery: predicting how compounds affect gene expression in different cell types without expensive experiments. Success could accelerate:
- Drug discovery and development
- Understanding of drug mechanisms
- Prediction of side effects
- Personalized medicine approaches

## Getting Started

1. Load the training data: `pd.read_parquet('de_train.parquet')`
2. Explore gene expression patterns across compounds and cell types
3. Build a model to predict test set responses
4. Format predictions using the sample submission template

The task requires predicting compound effects in B cells and Myeloid cells based on training data from multiple cell types and partial data from the target cell types.

