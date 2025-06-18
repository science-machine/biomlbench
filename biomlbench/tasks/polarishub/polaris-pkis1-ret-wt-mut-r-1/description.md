# pkis1-ret-wt-mut-r-1

## Overview

Multitask classification benchmark for RET wild type, mutant V804L, and mutant Y791F.

**Source:** [Polaris Hub - polaris/pkis1-ret-wt-mut-r-1](https://polarishub.io)  
**Task Type:** multi_task  
**Main Metric:** mean_squared_error

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `smiles`
- Target column: `{'RET_(V804L_mutant)', 'RET_(Y791F_mutant)', 'RET'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `polaris/pkis1-ret-wt-mut-r-1`.
Main metric: **mean_squared_error**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
