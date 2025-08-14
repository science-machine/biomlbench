# Cell-Cell Communication (Ligand-Target) Task

## Overview
Cell-cell Communication (ligand-target): Predicting cell-cell communication from single-cell transcriptomics data using supervised learning

## Task Description
The growing availability of single-cell data has sparked an increased interest in the inference of cell-cell communication (CCC), with an ever-growing number of computational tools developed for this purpose.

Different tools propose distinct preprocessing steps with diverse scoring functions, that are challenging to compare and evaluate. Furthermore, each tool typically comes with its own set of prior knowledge. To harmonize these, [Dimitrov et al, 2022](https://openproblems.bio/bibliography#dimitrov2022comparison) recently developed the [LIANA](https://github.com/saezlab/liana) framework, which was used as a foundation for this task.

The challenges in evaluating the tools are further exacerbated by the lack of a gold standard to benchmark the performance of CCC methods. In an attempt to address this, Dimitrov et al use alternative data modalities, including the spatial proximity of cell types and downstream cytokine activities, to generate an inferred ground truth. However, these modalities are only approximations of biological reality and come with their own assumptions and limitations. In time, the inclusion of more datasets with known ground truth interactions will become available, from which the limitations and advantages of the different CCC methods will be better understood.

**This subtask evaluates the methods' ability to predict interactions, the corresponding cytokines of which, are inferred to be active in the target cell types. This subtask focuses on the prediction of interactions from steady-state, or single-context, single-cell data.**

## Dataset: TNBC Data
We use the **Triple-Negative Breast Cancer (TNBC)** dataset from Wu et al., 2021. This dataset provides single-cell RNA-seq data from human breast cancer samples, offering insights into cell-cell communication within the tumor microenvironment.

### Dataset Overview:
- **Reference**: Wu et al., 2021 - Single-cell analysis of triple-negative breast cancer
- **Tissue**: Human breast cancer
- **Technology**: Single cell RNA sequencing (scRNA-seq)
- **Cell types**: Multiple tumor and immune cell types
- **Ground truth**: Inferred cytokine activities serving as proxy for cell-cell communication

## Data Format
All data is provided in H5AD (AnnData) format with the following structure:

### Input Data
- **obs**:
  - `label`: Cell type annotations
- **var**: Gene metadata
- **X**: Raw count matrix
- **uns**:
  - `ccc_train`: Training set DataFrame with columns:
    - `ligand`: Gene symbol of the ligand
    - `target`: Target cell type name
    - `response`: Binary (0/1) indicating interaction occurrence
  - `ccc_test_pairs`: Test set DataFrame with columns:
    - `ligand`: Gene symbol of the ligand
    - `target`: Target cell type name
    - (Note: response values are withheld for evaluation)
  - `target_organism`: NCBI taxonomy ID for species conversion

### Data Split
- The dataset is split 10/90 into training and test sets
- Training data includes ground truth responses for model development
- Test set responses are withheld and used for evaluation
- Stratified split ensures balanced class distribution

**NOTE**: We have purposely given a very small number of labels to you. This is because you are NOT expected to actually train a supervised model. The goal is to develop an unsupervised method / statistical approach to predict interactions. We have allowed you a small number of labels so that you can sanity check your method.

### Output Format
Methods must produce predictions for the test set interactions as a CSV file with:
- `ligand`: Gene symbol of the ligand (must match test set)
- `target`: Target cell type name (must match test set)
- `score`: Predicted interaction strength (-inf to +inf)

## Evaluation Metric
We evaluate methods using **Odds Ratio**, which:
- Compares true/false positive ratios in top-ranked predictions vs remaining interactions
- Uses top 5% of predictions by default
- Applies sigmoid transformation for normalization
- Quantifies association strength between method prioritization and positive interactions

Formula: OR = (TP × TN) / (FP × FN)
Where predictions are evaluated within top-ranked interactions.

## Input/Output Specification

### Input
Methods receive:
1. H5AD file (`tnbc_data.h5ad`) containing:
   - Cell type annotations and gene expression data
   - Training set with 50% of interactions (including ground truth responses)
   - Test set pairs (without responses) to predict
2. Prior knowledge ligand-receptor resource (`ligand_receptor_resource.csv.gz`):
   - OmniPath consensus database of known ligand-receptor interactions
   - Key columns: `source_genesymbol` (ligand), `target_genesymbol` (receptor)
   - Includes metadata like `secreted_intercell_source` (if ligand is secreted)
   - Aggregated from CellPhoneDB, CellChatDB, ICELLNET, connectomeDB2020, CellTalkDB
3. Gene symbol reference mapping (if available)

### Output
Methods must produce a CSV file with:
- Required columns: ligand, target, score
- Must contain exactly the same ligand-target pairs as in `ccc_test_pairs`
- Scores representing predicted interaction probability (0-1)

## Implementation Notes

### Using the Prior Knowledge Resource
The ligand-receptor resource (`ligand_receptor_resource.csv.gz`) provides biological constraints:
```python
import pandas as pd
lr_resource = pd.read_csv('ligand_receptor_resource.csv.gz')
# Key columns:
# - source_genesymbol: ligand gene
# - target_genesymbol: receptor gene  
# - secreted_intercell_source: True if ligand is secreted
```

It is the Omnipath database downloaded from LIANA's R package (Aug 2025).

## Important Considerations
1. **Prior knowledge integration**: Methods should leverage the ligand-receptor resource to:
   - Focus on biologically plausible interactions
   - Check if ligands are secreted (for cell-cell communication)
   - Handle multi-subunit complexes (subunits separated by `_`)
2. **Training data usage**: Use the provided training set to learn patterns
3. **Score interpretation**: Higher scores indicate stronger predicted interactions
4. **Cell type matching**: Target must be valid cell type from dataset
5. **Expression-based filtering**: Consider filtering interactions based on:
   - Ligand expression in any cell type
   - Receptor expression in the target cell type

