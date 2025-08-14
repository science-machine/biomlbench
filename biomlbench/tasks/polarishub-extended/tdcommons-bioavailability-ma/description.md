# bioavailability-ma

## Background
Oral bioavailability is defined as “the rate and extent to which the active ingredient or active moiety is absorbed from a drug product and becomes available at the site of action”.

## Description of readout
Task Description: Binary classification. Given a drug SMILES string, predict the activity of bioavailability.


## Data resource
**Reference**: [Ma, Chang-Ying, et al. “Prediction models of human plasma protein binding rate and oral bioavailability derived by using GA–CG–SVM method.” Journal of pharmaceutical and biomedical analysis 47.4-5 (2008): 677-682.](https://doi.org/10.1016/j.jpba.2008.03.023)


---

**Source:** [Polaris Hub - tdcommons/bioavailability-ma](https://polarishub.io)  
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

Uses official Polaris evaluation system with benchmark `tdcommons/bioavailability-ma`.
Main metric: **roc_auc**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
