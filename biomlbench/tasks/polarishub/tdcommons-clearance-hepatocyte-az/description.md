# clearance-hepatocyte-az

## Background
Drug clearance is defined as the volume of plasma cleared of a drug over a specified time period and it measures the rate at which the active drug is removed from the body. This is a dataset curated from ChEMBL database containing experimental results on intrinsic clearance, deposited from AstraZeneca. It contains clearance measures from two experiments types, hepatocyte and microsomes. As many studies [2] have shown various clearance outcomes given these two different types, we separate them.

## Description of readout
Task Description: Regression. Given a drug SMILES string, predict the activity of clearance.

## Data resource
**References**: [1] [Experimental in vitro Dmpk and physicochemical data on a set of publicly disclosed compounds](https://www.ebi.ac.uk/chembl/document_report_card/CHEMBL3301361/)

[2] [Mechanistic insights from comparing intrinsic clearance values between human liver microsomes and hepatocytes to guide drug design.](https://doi.org/10.1016/j.ejmech.2012.06.043)

---

**Source:** [Polaris Hub - tdcommons/clearance-hepatocyte-az](https://polarishub.io)  
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

Uses official Polaris evaluation system with benchmark `tdcommons/clearance-hepatocyte-az`.
Main metric: **spearmanr**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
