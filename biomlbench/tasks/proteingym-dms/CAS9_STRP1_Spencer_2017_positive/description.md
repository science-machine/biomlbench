
# Proteingym-DMS dataset: CAS9_STRP1_Spencer_2017_positive

## Description

This dataset is part of the ProteinGym DMS benchmark, which contains deep mutational scanning datasets that measure
protein fitness (in different contexts) for sequence variants of a wide range of proteins. This dataset contains
only single-substitution variants for the protein Streptococcus pyogenes Cas9 from the organism Streptococcus pyogenes. This protein has Uniprot ID: CAS9_STRP1. 

The DMS selection assay was described as follows: 

    Count of mutation where survival depends on expression of cas9 and correct cleavage

It was categorised as measuring the following (general) fitness attribute: Activity. Higher scores indicate better fitness.

The source publication for this dataset is titled: 

"Deep mutational scanning of S. pyogenes Cas9 reveals important functional domains"

and can be accessed at the following DOI: 10.1038/s41598-017-17081-y.

## Objective

The objective of this benchmark is to train a model that can predict the fitness of unseen single-substitution sequence variants of Streptococcus pyogenes Cas9.
These variants may be in regions of the protein where mutations were not present in the training set.

## Data Format

This task uses CSV files to store the amino acid sequence data and fitness scores.
- ID column: `id` (this is just an index column)
- Sequence column: `sequence`
- Fitness column: `fitness_score`

The first row of the CSV file contains the wild-type sequence in the `sequence` field and a missing `fitness_score` field.

## Files

- `train.csv`: Training data with sequences and fitness scores
- `test_features.csv`: Test features with ID column and sequence column - these fitness scores must be predicted by your model
- `sample_submission.csv`: Example submission format with ID column and fitness score column

## Evaluation

Your model will be evaluated on the Spearman correlation between the predicted fitness scores and the true fitness scores.
