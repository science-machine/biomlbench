# rna-expression-prediction-dataset-task



---

**Source:** [Polaris Hub - vishrut64/rna-expression-prediction-dataset-task](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** mean_absolute_error

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `cleaned_sequence`
- Target column: `{'expression'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `vishrut64/rna-expression-prediction-dataset-task`.
Main metric: **mean_absolute_error**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
