# adme-fang-reg-v1

## Overview

A multitask benchmark for six ADME endpoints, utilizing a shared random split. 

**Source:** [Polaris Hub - biogen/adme-fang-reg-v1](https://polarishub.io)  
**Task Type:** multi_task  
**Main Metric:** mean_squared_error

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `MOL_smiles`
- Target column: `{'LOG_RLM_CLint', 'LOG_SOLUBILITY', 'LOG_HPPB', 'LOG_MDR1-MDCK_ER', 'LOG_HLM_CLint', 'LOG_RPPB'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `biogen/adme-fang-reg-v1`.
Main metric: **mean_squared_error**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
