# molprop-leadlike-250k-reg-v2

## Overview

A multitask benchmark designed to predict nine molecular properties for 250,000 compounds sourced from ZINC15, with a focus on molecular representation.

**Source:** [Polaris Hub - polaris/molprop-leadlike-250k-reg-v2](https://polarishub.io)  
**Task Type:** multi_task  
**Main Metric:** mean_squared_error

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `MOL_smiles`
- Target column: `{'mw', 'n_charged_atoms', 'fsp3', 'n_aromatic_rings', 'tpsa', 'clogp', 'n_rotatable_bonds', 'refractivity', 'formal_charge'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `polaris/molprop-leadlike-250k-reg-v2`.
Main metric: **mean_squared_error**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
