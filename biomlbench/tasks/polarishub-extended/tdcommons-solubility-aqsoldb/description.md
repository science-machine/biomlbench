# solubility-aqsoldb

## Background
Aqeuous solubility measures a drug's ability to dissolve in water. Poor water solubility could lead to slow drug absorptions, inadequate bioavailablity and even induce toxicity. More than 40% of new chemical entities are not soluble.

## Description of readout
Task Description: Regression. Given a drug SMILES string, predict the activity of solubility.

## Data resource
**Reference**: [1] [AqSolDB, a curated reference set of aqueous solubility and 2D descriptors for a diverse set of compounds](https://www.nature.com/articles/s41597-019-0151-1)

---

**Source:** [Polaris Hub - tdcommons/solubility-aqsoldb](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** mean_absolute_error

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

Uses official Polaris evaluation system with benchmark `tdcommons/solubility-aqsoldb`.
Main metric: **mean_absolute_error**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
