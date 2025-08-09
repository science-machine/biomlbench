# bbb-martins

## Background
As a membrane separating circulating blood and brain extracellular fluid, the blood-brain barrier (BBB) is the protection layer that blocks most foreign drugs. Thus the ability of a drug to penetrate the barrier to deliver to the site of action forms a crucial challenge in development of drugs for central nervous system From MoleculeNet.

## Description of readout
Task Description: Binary classification. Given a drug SMILES string, predict the activity of BBB.

## Data resource
**Reference**: [1] [A Bayesian approach to in silico blood-brain barrier penetration modeling.](https://pubs.acs.org/doi/10.1021/ci300124c)

[2] [MoleculeNet: a benchmark for molecular machine learning.](https://pubs.rsc.org/en/content/articlelanding/2018/sc/c7sc02664a)

---

**Source:** [Polaris Hub - tdcommons/bbb-martins](https://polarishub.io)  
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

Uses official Polaris evaluation system with benchmark `tdcommons/bbb-martins`.
Main metric: **roc_auc**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
