# Predict Modality Task

## Overview
Predict Modality: Predicting the profiles of one modality (e.g. protein abundance) from another (e.g. mRNA expression).

## Task Description
Experimental techniques to measure multiple modalities within the same single cell are increasingly becoming available. 
The demand for these measurements is driven by the promise to provide a deeper insight into the state of a cell. 
Yet, the modalities are also intrinsically linked. We know that DNA must be accessible (ATAC data) to produce mRNA 
(expression data), and mRNA in turn is used as a template to produce protein (protein abundance). These processes 
are regulated often by the same molecules that they produce: for example, a protein may bind DNA to prevent the production 
of more mRNA. Understanding these regulatory processes would be transformative for synthetic biology and drug target discovery. 
Any method that can predict a modality from another must have accounted for these regulatory processes, but the demand for 
multi-modal data shows that this is not trivial.

## Dataset: BMMC CITE-seq
This task uses a Bone Marrow Mononuclear Cells (BMMC) CITE-seq dataset. CITE-seq (Cellular Indexing of Transcriptomes and Epitopes by Sequencing) is a method that simultaneously measures:
- **Modality 1 (RNA)**: Gene expression profiles (mRNA abundance)
- **Modality 2 (Protein)**: Surface protein abundance measured using oligonucleotide-labeled antibodies

The dataset contains:
- **Training data**: Paired RNA and protein measurements from the same cells
- **Test data**: RNA measurements only (participants must predict the corresponding protein measurements)

Key statistics:
- Number of genes (RNA): ~13,000
- Number of proteins: ~134
- Number of cells in training set: varies
- Number of cells in test set: varies

## Data Format
All data is provided in H5AD (AnnData) format:
- `train_mod1.h5ad`: Training RNA expression data
- `train_mod2.h5ad`: Training protein abundance data (matched cells with train_mod1)
- `test_mod1.h5ad`: Test RNA expression data

The data has been preprocessed:
- Log-transformed
- Normalized to counts per 10,000 (CP10K)
- Expression values are stored in the `layers['normalized']` attribute (not in `.X`)

## Evaluation Metric
The task uses **Root Mean Squared Error (RMSE)** calculated across all protein-cell pairs. The score is computed as:
```
RMSE = sqrt(mean((y_true - y_pred)^2))
```
where:
- `y_true`: Ground truth protein abundance values
- `y_pred`: Predicted protein abundance values

Lower RMSE indicates better performance.

## Input/Output Specification

### Input
Participants receive:
1. `train_mod1.h5ad`: RNA expression data for training cells
2. `train_mod2.h5ad`: Protein abundance data for the same training cells
3. `test_mod1.h5ad`: RNA expression data for test cells

### Output
Participants must submit a CSV file (`submission.csv`) with the following columns:
- `cell_id`: Cell identifier matching the order in test_mod1
- Remaining columns: One column per protein, containing predicted abundance values

The column names must exactly match the protein names (var_names) from train_mod2.

### Sample Submission Format
```csv
cell_id,CD3,CD4,CD8,CD19,CD20,...
AAACCCAAGAAACTAG-1,0.123,0.456,0.789,0.012,0.345,...
AAACCCAAGAAACTAT-1,0.234,0.567,0.890,0.123,0.456,...
...
```

