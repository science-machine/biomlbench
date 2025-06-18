# cyp2c9-veith

## Background
The CYP P450 genes are involved in the formation and breakdown (metabolism) of various molecules and chemicals within cells. Specifically, the CYP P450 2C9 plays a major role in the oxidation of both xenobiotic and endogenous compounds.

## Description of readout
Task Description: Binary Classification. Given a drug SMILES string, predict CYP2C9 inhibition.

## Data resource
**References**: [1] [Comprehensive Characterization of Cytochrome P450 Isozyme Selectivity across Chemical Libraries
](https://www.nature.com/articles/nbt.1581)

---

**Source:** [Polaris Hub - tdcommons/cyp2c9-veith](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** pr_auc

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

Uses official Polaris evaluation system with benchmark `tdcommons/cyp2c9-veith`.
Main metric: **pr_auc**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
