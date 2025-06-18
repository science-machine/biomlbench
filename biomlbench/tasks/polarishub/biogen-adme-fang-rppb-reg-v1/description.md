# adme-fang-RPPB-reg-v1

## Overview

Single task benchmark for ADME property LOG_RPPB

**Source:** [Polaris Hub - biogen/adme-fang-rppb-reg-v1](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** pearsonr

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `MOL_smiles`
- Target column: `{'LOG_RPPB'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `biogen/adme-fang-rppb-reg-v1`.
Main metric: **pearsonr**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
