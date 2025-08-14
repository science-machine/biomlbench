# hia-hou

## Background
When a drug is orally administered, it needs to be absorbed from the human gastrointestinal system into the bloodstream of the human body. This ability of absorption is called human intestinal absorption (HIA) and it is crucial for a drug to be delivered to the target.

## Description of readout
Task Description: Binary classification. Given a drug SMILES string, predict the activity of HIA.


## Data resource
**Reference**: [Hou T et al. ADME evaluation in drug discovery. 7. Prediction of oral absorption by correlation and classification. J Chem Inf Model. 2007;47(1):208-218.](https://pubmed.ncbi.nlm.nih.gov/17238266/)


---

**Source:** [Polaris Hub - tdcommons/hia-hou](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** roc_auc

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

Uses official Polaris evaluation system with benchmark `tdcommons/hia-hou`.
Main metric: **roc_auc**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
