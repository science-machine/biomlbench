# pkis2-egfr-wt-c-1

## Overview

Single task classification benchmark for kinase EGFR wild type.

**Source:** [Polaris Hub - polaris/pkis2-egfr-wt-c-1](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** pr_auc

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `smiles`
- Target column: `{'CLASS_EGFR'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `polaris/pkis2-egfr-wt-c-1`.
Main metric: **pr_auc**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
