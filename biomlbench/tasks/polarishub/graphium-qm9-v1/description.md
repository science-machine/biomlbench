# qm9-v1

## Overview

A multitask regression benchmark for QM9 dataset

**Source:** [Polaris Hub - graphium/qm9-v1](https://polarishub.io)  
**Task Type:** multi_task  
**Main Metric:** mean_squared_error

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `smiles`
- Target column: `{'g298_atom', 'gap', 'r2', 'u298_atom', 'homo', 'h298_atom', 'u0', 'C', 'B', 'zpve', 'h298', 'alpha', 'mu', 'g298', 'lumo', 'u0_atom', 'A', 'u298', 'cv'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `graphium/qm9-v1`.
Main metric: **mean_squared_error**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
