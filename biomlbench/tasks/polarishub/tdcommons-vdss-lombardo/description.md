# vdss-lombardo

## Background
The volume of distribution at steady state (VDss) measures the degree of a drug's concentration in body tissue compared to concentration in blood. Higher VD indicates a higher distribution in the tissue and usually indicates the drug with high lipid solubility, low plasma protein binidng rate.

## Description of readout
Task Description: Regression. Given a drug SMILES string, predict the volume of distributon.

## Data resource
**Reference**: [1] [In Silico Prediction of Volume of Distribution in Humans. Extensive Data Set and the Exploration of Linear and Nonlinear Methods Coupled with Molecular Interaction Fields Descriptors](https://pubs.acs.org/doi/10.1021/acs.jcim.6b00044)

---

**Source:** [Polaris Hub - tdcommons/vdss-lombardo](https://polarishub.io)  
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

Uses official Polaris evaluation system with benchmark `tdcommons/vdss-lombardo`.
Main metric: **spearmanr**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
