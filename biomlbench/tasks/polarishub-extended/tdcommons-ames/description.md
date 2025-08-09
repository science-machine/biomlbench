# ames

## Background
Mutagenicity means the ability of a drug to induce genetic alterations. Drugs that can cause damage to the DNA can result in cell death or other severe adverse effects. Nowadays, the most widely used assay for testing the mutagenicity of compounds is the Ames experiment which was invented by a professor named Ames. The Ames test is a short-term bacterial reverse mutation assay detecting a large number of compounds which can induce genetic damage and frameshift mutations. The dataset is aggregated from four papers.

## Description of readout
Task Description: Binary classification. Given a drug SMILES string, predict whether it is mutagenic (1) or not mutagenic (0).

## Data resource
**Reference**: [1] [In silico Prediction of Chemical Ames Mutagenicity](https://pubs.acs.org/doi/10.1021/ci300400a)

---

**Source:** [Polaris Hub - tdcommons/ames](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** roc_auc

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `Drug`
- Target column: `{'Y'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `tdcommons/ames`.
Main metric: **roc_auc**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
