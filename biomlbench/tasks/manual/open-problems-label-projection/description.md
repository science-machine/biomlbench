# Label Projection Task

## Overview
Label projection: Automated cell type annotation from rich, labeled reference data

## Task Description
A major challenge for integrating single cell datasets is creating matching
cell type annotations for each cell. One of the most common strategies for
annotating cell types is referred to as
["cluster-then-annotate"](https://www.nature.com/articles/s41576-018-0088-9)
whereby cells are aggregated into clusters based on feature similarity and
then manually characterized based on differential gene expression or previously
identified marker genes. Recently, methods have emerged to build on this
strategy and annotate cells using
[known marker genes](https://www.nature.com/articles/s41592-019-0535-3).
However, these strategies pose a difficulty for integrating atlas-scale
datasets as the particular annotations may not match.

To ensure that the cell type labels in newly generated datasets match
existing reference datasets, some methods align cells to a previously
annotated [reference dataset](https://academic.oup.com/bioinformatics/article/35/22/4688/54802990)
and then _project_ labels from the reference to the new dataset.

Here, we compare methods for annotation based on a reference dataset.
The datasets consist of two or more samples of single cell profiles that
have been manually annotated with matching labels. These datasets are then
split into training and test batches, and the task of each method is to
train a cell type classifer on the training set and project those labels
onto the test set.

## Dataset: Diabetic Kidney Disease (DKD)
We use the **Diabetic Kidney Disease** dataset from the CellxGene Census, which provides comprehensive single-cell data from human kidney cortex samples. This dataset includes cells from donors with and without diabetic kidney disease, offering insights into disease-associated cellular changes.

### Dataset Overview:
- **Reference**: Wilson et al., 2022 - "Multimodal single cell sequencing implicates chromatin accessibility and genetic background in diabetic kidney disease progression"
- **Tissue**: Human kidney cortex
- **Condition**: Samples from donors with and without diabetic kidney disease
- **Technology**: Single nucleus RNA sequencing (snRNA-seq)
- **Cell types**: Multiple kidney cell types including:
  - Proximal tubule cells (including injured PT_VCAM1+ cells)
  - Distal tubule cells
  - Collecting duct cells
  - Glomerular cells (podocytes, endothelial, mesangial)
  - Immune cells
  - Stromal cells

## Data Format
All data is provided in H5AD (AnnData) format with the following structure:

### Training Data (`train.h5ad`)
- **obs**: 
  - `label`: Cell type annotations
  - `batch`: Batch/sample information
- **var**: Gene metadata (feature_id, feature_name, hvg, hvg_score)
- **layers**:
  - `counts`: Raw count data
  - `normalized`: Log-normalized expression (log(CP10k+1))
- **obsm**: 
  - `X_pca`: Pre-computed PCA embeddings

### Test Data (`test.h5ad`)
- Same structure as training data but without `label` column
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
s3://openproblems-data/resources/task_label_projection/datasets/cellxgene_census/dkd/log_cp10k/
``` 