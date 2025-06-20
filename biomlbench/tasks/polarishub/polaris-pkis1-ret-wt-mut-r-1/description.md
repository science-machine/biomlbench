# pkis1-ret-wt-mut-r-1

## Background
**RET (Rearranged during Transfection)** is a proto-oncogene that codes for a receptor tyrosine kinase. This means it produces a protein that plays a role in signaling pathways within cells, particularly related to cell growth and differentiation. When activated, RET helps regulate cell survival, proliferation, and differentiation. Mutations or alterations in the RET gene can lead to uncontrolled cell growth and potentially the development of cancer.

## Benchmarking
**RET wild type**: In some cases, targeting both mutant and wild-type RET together can be more effective than targeting only one form as Combination Therapies. In certain cancer types, such as some subtypes of non-small cell lung cancer (NSCLC), the RET signaling pathway can interact with other oncogenic pathways, such as the EGFR (epidermal growth factor receptor) pathway. Targeting both pathways simultaneously might offer a synergistic effect and improve treatment outcomes. 

**RET-V804L**: The V804L mutation causes a structural change in the RET protein, resulting in its continuous activation, even in the absence of ligand binding. The mutation leads to uncontrolled cell growth and division, contributing to oncogenesis. The V804L mutation in RET has been identified in various cancer types, particularly in thyroid cancers, including papillary thyroid carcinoma (PTC) and medullary thyroid carcinoma (MTC). It is often associated with aggressive tumor behavior and resistance to conventional therapies. 

**RET-Y791F**: The Y791F mutation disrupts a crucial phosphorylation site within the RET protein. Consequently, the tyrosine at position 791 cannot be phosphorylated effectively or at all. This disruption interferes with the normal signaling pathways that rely on this particular phosphorylation event, leading to dysregulated downstream signaling, such as MAPK and PI3K. The Y791F mutation, by abrogating this phosphorylation site, alters these cellular responses and may influence tumor development and progression. The Y791F mutation is found in the RET protein, which is implicated in several types of cancer, particularly medullary thyroid carcinoma (MTC).

**The goal** of this benchmark is to select the best predictive model for 
- Selectivity towards the mutants
- Optimization of the bioactivity % inhibition.
- Discovery of potential hits in new chemical space.


## Description of readout 
**Description of readout:**
- **Readouts**: `RET`, `RET_(V804L_mutant)`, `RET_(Y791F_mutant)`
- **Bioassay readout**: Percentage of inhibition.
- **Optimization objective**: Higher value (higher %inhibition ).
- **Number of data points**: train:  279 test:  85

## Data resource: 
- **Reference**: [PKIS1](https://pubmed.ncbi.nlm.nih.gov/26501955)


## Train/test split
Given the benchmarking goal, a scaffold-based splitting approach was applied to ensure training and test sets contain distinct chemical structures while maintaining the diversity of scaffolds.

**Distribution of the train/test in the chemical space**
![image](https://storage.googleapis.com/polaris-public/datasets/kinases/ret/figures/drewry_ret_wt_v804l_y791f_v1_tnse_scaffold_split.png)

**For more details of this benchmark** -> [notebook](https://github.com/polaris-hub/polaris-recipes/blob/main/03_Kinases/RET)



---

**Source:** [Polaris Hub - polaris/pkis1-ret-wt-mut-r-1](https://polarishub.io)  
**Task Type:** multi_task  
**Main Metric:** mean_squared_error

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `smiles`
- Target column: `{'RET_(V804L_mutant)', 'RET', 'RET_(Y791F_mutant)'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `polaris/pkis1-ret-wt-mut-r-1`.
Main metric: **mean_squared_error**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
