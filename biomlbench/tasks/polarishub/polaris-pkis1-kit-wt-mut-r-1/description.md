# pkis1-kit-wt-mut-r-1

## Overview

A multitask regression benchmark for KIT wild type, T670I mutant and KV560G_mutant.

**Source:** [Polaris Hub - polaris/pkis1-kit-wt-mut-r-1](https://polarishub.io)  
**Task Type:** multi_task  
**Main Metric:** mean_squared_error

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `smiles`
- Target column: `{'KIT_(V560G_mutant)', 'KIT_(T6701_mutant)', 'KIT'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `polaris/pkis1-kit-wt-mut-r-1`.
Main metric: **mean_squared_error**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
