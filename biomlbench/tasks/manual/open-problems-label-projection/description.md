# Label Projection Task

## Overview
Label projection is a fundamental task in single-cell genomics where the goal is to automatically annotate cell types in new datasets based on existing reference atlases. This task evaluates methods for transferring cell type labels from a well-annotated reference dataset to an unlabeled query dataset.

## Task Description
Given a training dataset with cell type annotations and an unlabeled test dataset, methods must predict the cell type for each cell in the test set. This is a multi-class classification problem that requires:

1. Learning cell type signatures from the reference data
2. Handling batch effects between datasets
3. Projecting labels accurately to new cells

## Dataset: Mouse Pancreas Atlas
We use the **Mouse Pancreas Atlas** from the CellxGene Census, which provides a comprehensive single-cell atlas of mouse pancreatic cells. This dataset includes:

- **Training set**: 292,354 cells with known cell type labels
- **Test set**: 9,442 cells requiring label prediction
- **Cell types**: Multiple pancreatic cell types including:
  - Type A, B, D pancreatic cells (alpha, beta, delta cells)
  - Pancreatic PP cells
  - Pancreatic ductal cells
  - Pancreatic stellate cells
  - Endothelial cells
  - Hematopoietic cells
  - And others

## Data Format
All data is provided in H5AD (AnnData) format with the following structure:

### Training Data (`train.h5ad`)
- **Dimensions**: 292,354 cells × 24,923 genes
- **obs**: 
  - `label`: Cell type annotations
  - `batch`: Batch/study information
- **var**: Gene metadata (feature_id, feature_name, hvg, hvg_score)
- **layers**:
  - `counts`: Raw count data
  - `normalized`: Log-normalized expression (log(CP10k+1))
- **obsm**: 
  - `X_pca`: Pre-computed PCA embeddings

### Test Data (`test.h5ad`)
- **Dimensions**: 9,442 cells × 24,923 genes
- **Structure**: Same as training data but without `label` column
- Methods must predict labels for these cells

## Evaluation Metric
We evaluate methods using **F1-weighted score**, which:
- Calculates F1 score for each cell type class
- Weights by the number of true instances for each class
- Provides a balanced measure that accounts for class imbalance

Formula: F1-weighted = Σ(n_i / N) × F1_i

Where:
- n_i = number of samples in class i
- N = total number of samples
- F1_i = F1 score for class i

## Input/Output Specification

### Input
Methods receive:
1. `train.h5ad`: Training data with cell type labels
2. `test.h5ad`: Test data without labels

### Output
Methods must produce a prediction file with:
- Cell identifiers matching the test set
- Predicted cell type label for each cell
- Labels must exactly match the vocabulary from the training set

## Implementation Notes
- Methods can use any layer (`counts` or `normalized`) or the pre-computed PCA
- Batch information is provided but handling batch effects is up to each method
- Gene sets are pre-aligned between train and test sets
- All cells must receive a prediction (no missing values allowed)

## Baseline Performance
The task includes control methods for comparison:
- **Random Labels**: Random assignment from training label distribution
- **Majority Vote**: Assigns the most frequent cell type to all cells
- **True Labels**: Oracle performance (upper bound)

## Important Considerations
1. **Exact label matching**: Predicted labels must exactly match training set vocabulary
2. **No missing predictions**: Every test cell must have a prediction
3. **Deterministic results**: Methods should produce consistent results across runs
4. **Resource constraints**: Methods should complete within reasonable time/memory limits

## Data Access
Data is hosted on the OpenProblems S3 bucket:
```
s3://openproblems-data/resources/task_label_projection/datasets/cellxgene_census/mouse_pancreas_atlas/log_cp10k/
``` 