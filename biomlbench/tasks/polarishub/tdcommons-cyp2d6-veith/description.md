# cyp2d6-veith

## Background
The CYP P450 genes are involved in the formation and breakdown (metabolism) of various molecules and chemicals within cells. Specifically, CYP2D6 is primarily expressed in the liver. It is also highly expressed in areas of the central nervous system, including the substantia nigra.

## Description of readout
Task Description: Binary Classification. Given a drug SMILES string, predict CYP2D6 inhibition.

## Data resource
**References**: [1] [Comprehensive Characterization of Cytochrome P450 Isozyme Selectivity across Chemical Libraries
](https://www.nature.com/articles/nbt.1581)

---

**Source:** [Polaris Hub - tdcommons/cyp2d6-veith](https://polarishub.io)  
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

Uses official Polaris evaluation system with benchmark `tdcommons/cyp2d6-veith`.
Main metric: **pr_auc**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
