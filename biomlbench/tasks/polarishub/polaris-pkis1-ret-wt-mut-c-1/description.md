# pkis1-ret-wt-mut-c-1

## Overview

Multitask classification benchmark for RET wild type, mutant V804L, and mutant Y791F.

**Source:** [Polaris Hub - polaris/pkis1-ret-wt-mut-c-1](https://polarishub.io)  
**Task Type:** multi_task  
**Main Metric:** pr_auc

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `smiles`
- Target column: `{'CLASS_RET_(V804L_mutant)', 'CLASS_RET', 'CLASS_RET_(Y791F_mutant)'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `polaris/pkis1-ret-wt-mut-c-1`.
Main metric: **pr_auc**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
