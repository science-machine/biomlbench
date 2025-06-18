# clearance-microsome-az

## Overview



**Source:** [Polaris Hub - tdcommons/clearance-microsome-az](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** spearmanr

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `Drug`
- Target column: `{'Y'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `tdcommons/clearance-microsome-az`.
Main metric: **spearmanr**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
