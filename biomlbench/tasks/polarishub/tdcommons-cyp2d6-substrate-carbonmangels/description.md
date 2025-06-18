# cyp2d6-substrate-carbonmangels

## Background
CYP2D6 is primarily expressed in the liver. It is also highly expressed in areas of the central nervous system, including the substantia nigra. TDC used a dataset from [1], which merged information on substrates and nonsubstrates from six publications.

## Description of readout
Task Description: Binary Classification. Given a drug SMILES string, predict if it is a substrate to the enzyme.

## Data resource
**References**: [1] [Selecting relevant descriptors for classification by bayesian estimates: a comparison with decision trees and support vector machines approaches for disparate data sets.](https://doi.org/10.1002/minf.201100069)

[2] [admetSAR: a comprehensive source and free tool for assessment of chemical ADMET properties.](https://pubs.acs.org/doi/10.1021/ci300367a)

---

**Source:** [Polaris Hub - tdcommons/cyp2d6-substrate-carbonmangels](https://polarishub.io)  
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

Uses official Polaris evaluation system with benchmark `tdcommons/cyp2d6-substrate-carbonmangels`.
Main metric: **pr_auc**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
