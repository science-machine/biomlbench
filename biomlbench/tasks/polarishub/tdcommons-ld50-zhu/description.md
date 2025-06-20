# ld50-zhu

## Background
Acute toxicity LD50 measures the most conservative dose that can lead to lethal adverse effects. The higher the dose, the more lethal of a drug. This dataset is kindly provided by the authors of [1]

## Description of readout
Task Description: Regression. Given a drug SMILES string, predict its acute toxicity.

## Data resource
**Reference**: [1] [Quantitative Structure−Activity Relationship Modeling of Rat Acute Toxicity by Oral Exposure](https://pubs.acs.org/doi/10.1021/tx900189p)


---

**Source:** [Polaris Hub - tdcommons/ld50-zhu](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** mean_absolute_error

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

Uses official Polaris evaluation system with benchmark `tdcommons/ld50-zhu`.
Main metric: **mean_absolute_error**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
