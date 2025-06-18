# posebusters_v1

## Overview

Docking task benchmark for 428 proteins and ligands.

**Source:** [Polaris Hub - polaris/posebusters-v1](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** rmsd_coverage

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `protein`
- Target column: `{'ligand'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `polaris/posebusters-v1`.
Main metric: **rmsd_coverage**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
