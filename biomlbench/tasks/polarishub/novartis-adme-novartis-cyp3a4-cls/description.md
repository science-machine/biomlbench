# adme-novartis-cyp3a4-cls

## Overview

Single classification task benchmark for CYP3A4 log_kobs

**Source:** [Polaris Hub - novartis/adme-novartis-cyp3a4-cls](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** balanced_accuracy

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `MOL_smiles`
- Target column: `{'CLS_log_kobs'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `novartis/adme-novartis-cyp3a4-cls`.
Main metric: **balanced_accuracy**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
