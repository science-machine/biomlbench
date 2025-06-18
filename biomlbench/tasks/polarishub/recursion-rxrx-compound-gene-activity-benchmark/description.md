# RxRx-compound-gene-activity-benchmark

## Overview

This benchmark evaluates the zero-shot prediction of compound-gene activity. It is originally designed to be used with the rxrx.ai dataset releases in phenomics.

**Source:** [Polaris Hub - recursion/rxrx-compound-gene-activity-benchmark](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** pr_auc{"aggregation": "mean", "default": null, "group_by": "treatment", "on_error": "ignore"}

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `gene_symbol`
- Target column: `{'active'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `recursion/rxrx-compound-gene-activity-benchmark`.
Main metric: **pr_auc{"aggregation": "mean", "default": null, "group_by": "treatment", "on_error": "ignore"}**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
