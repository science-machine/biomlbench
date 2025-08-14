# ppbr-az

## Background
The human plasma protein binding rate (PPBR) is expressed as the percentage of a drug bound to plasma proteins in the blood. This rate strongly affect a drug's efficiency of delivery. The less bound a drug is, the more efficiently it can traverse and diffuse to the site of actions. From a ChEMBL assay deposited by AstraZeneca.

## Description of readout
Task Description: Regression. Given a drug SMILES string, predict the rate of PPBR.

## Data resource
**Reference**: [1] [Experimental in vitro Dmpk and physicochemical data on a set of publicly disclosed compounds](https://www.ebi.ac.uk/chembl/document_report_card/CHEMBL3301361/)

---

**Source:** [Polaris Hub - tdcommons/ppbr-az](https://polarishub.io)  
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

Uses official Polaris evaluation system with benchmark `tdcommons/ppbr-az`.
Main metric: **mean_absolute_error**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
