# moleculeace-CHEMBL1862-Ki

## Background
This dataset comprises collected and curated bioactivity data for the target Tyrosine protein kinase Abl familyfrom ChEMBL assay CHEMBL1862, utilized to evaluate the performance of various machine learning algorithms on activity cliffs. We employed classical machine learning methods combined with common molecular descriptors, as well as neural networks based on unstructured molecular data such as molecular graphs or SMILES strings.

Activity cliffs are molecules with small differences in structure but large differences in potency. Activity cliffs play an important role in drug discovery, but the bioactivity of activity cliff compounds are notoriously difficult to predict.

## Description of readouts
- `exp_mean [nM]`: Inhibition [Inhibitory Constant, Ki]
- `y`: Negative of log transform of the bioactivity value.
- `split`: Train-test split based on activity cliff.

## Data resource:
- [Exposing the Limitations of Molecular Machine Learning with Activity Cliffs
](https://pubs.acs.org/doi/10.1021/acs.jcim.2c01073)
- Github: https://github.com/molML/MoleculeACE

---

**Source:** [Polaris Hub - molecularml/moleculeace-chembl1862-ki](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** mean_squared_error

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `smiles`
- Target column: `{'y'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `molecularml/moleculeace-chembl1862-ki`.
Main metric: **mean_squared_error**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
