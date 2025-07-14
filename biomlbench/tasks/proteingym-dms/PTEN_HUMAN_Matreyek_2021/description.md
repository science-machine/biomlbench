
# Proteingym-DMS dataset: PTEN_HUMAN_Matreyek_2021

## Description

This dataset is part of the ProteinGym DMS benchmark, which contains deep mutational scanning datasets that measure
protein fitness (in different contexts) for sequence variants of a wide range of proteins. This dataset contains
only single-substitution variants for the protein PTEN from the organism Homo sapiens. This protein has Uniprot ID: PTEN_HUMAN. 

The DMS selection assay was described as follows: 

    Protein abundance (facs sorting for abundance of gfp-fused target)

It was categorised as measuring the following (general) fitness attribute: Expression. Higher scores indicate better fitness.

The source publication for this dataset is titled: 

"Integrating thousands of PTEN variant activity and abundance measurements reveals variant subgroups and new dominant negatives in cancers"

and can be accessed at the following DOI: 10.1186/s13073-021-00984-x.

## Objective

The objective of this benchmark is to train a model that can predict the fitness of unseen single-substitution sequence variants of PTEN.
These variants may be in regions of the protein where mutations were not present in the training set.

## Data Format

This task uses CSV files to store the amino acid sequence data and fitness scores.
- ID column: `id` (this is just an index column)
- Sequence column: `sequence`
- Fitness column: `fitness_score`

## Files

- `train.csv`: Training data with sequences and fitness scores
- `test_features.csv`: Test features with ID column and sequence column - these fitness scores must be predicted by your model
- `sample_submission.csv`: Example submission format with ID column and fitness score column

## Evaluation

Your model will be evaluated on the Spearman correlation between the predicted fitness scores and the true fitness scores.
