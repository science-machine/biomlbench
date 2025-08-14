# pgp-broccatelli

## Background
P-glycoprotein (Pgp) is an ABC transporter protein involved in intestinal absorption, drug metabolism, and brain penetration, and its inhibition can seriously alter a drug's bioavailability and safety. In addition, inhibitors of Pgp can be used to overcome multidrug resistance.

## Description of readout
Task Description: Binary classification. Given a drug SMILES string, predict the activity of Pgp inhibition.

## Data resource
**Reference**: [Broccatelli et al., A Novel Approach for Predicting P-Glycoprotein (ABCB1) Inhibition Using Molecular Interaction Fields. Journal of Medicinal Chemistry, 2011 54 (6), 1740-1751](https://pubs.acs.org/doi/abs/10.1021/jm101421d)


---

**Source:** [Polaris Hub - tdcommons/pgp-broccatelli](https://polarishub.io)  
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

Uses official Polaris evaluation system with benchmark `tdcommons/pgp-broccatelli`.
Main metric: **roc_auc**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
