# lipophilicity-astrazeneca

## Background
Lipophilicity measures the ability of a drug to dissolve in a lipid (e.g. fats, oils) environment. High lipophilicity often leads to high rate of metabolism, poor solubility, high turn-over, and low absorption. From MoleculeNet.

## Description of readout
Task Description: Regression. Given a drug SMILES string, predict the activity of lipophilicity.


## Data resource
**Reference**: 

- [AstraZeneca. Experimental in vitro Dmpk and physicochemical data on a set of publicly disclosed compounds (2016)](https://doi.org/10.6019/chembl3301361)

- [Wu, Zhenqin, et al. “MoleculeNet: a benchmark for molecular machine learning.” Chemical science 9.2 (2018): 513-530.](https://pubs.rsc.org/--/content/articlehtml/2018/sc/c7sc02664a)


---

**Source:** [Polaris Hub - tdcommons/lipophilicity-astrazeneca](https://polarishub.io)  
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

Uses official Polaris evaluation system with benchmark `tdcommons/lipophilicity-astrazeneca`.
Main metric: **mean_absolute_error**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
