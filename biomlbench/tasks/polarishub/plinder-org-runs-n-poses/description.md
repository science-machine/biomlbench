# Plinder-Org - Runs-N-Poses

## Overview

Benchmark from Polaris Hub: plinder-org/runs-n-poses

**Source:** [Polaris Hub - plinder-org/runs-n-poses](https://polarishub.io)  
**Task Type:** unknown  
**Main Metric:** unknown

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `smiles`
- Target column: `target` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `plinder-org/runs-n-poses`.
Main metric: **unknown**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
