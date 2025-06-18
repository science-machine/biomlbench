# half-life-obach

## Background
Half life of a drug is the duration for the concentration of the drug in the body to be reduced by half. It measures the duration of actions of a drug. This dataset is from [1] and we obtain the deposited version under CHEMBL assay 1614674.

## Description of readout
Task Description: Regression. Given a drug SMILES string, predict the half life duration.

## Data resource
**References**: [1] [Trend Analysis of a Database of Intravenous Pharmacokinetic Parameters in Humans for 670 Drug Compounds
](https://doi.org/10.1124/dmd.108.020479)

---

**Source:** [Polaris Hub - tdcommons/half-life-obach](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** spearmanr

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

Uses official Polaris evaluation system with benchmark `tdcommons/half-life-obach`.
Main metric: **spearmanr**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
