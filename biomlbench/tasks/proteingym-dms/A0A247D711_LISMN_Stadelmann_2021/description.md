
# Proteingym-DMS dataset: A0A247D711_LISMN_Stadelmann_2021

## Description

This dataset is part of the ProteinGym DMS benchmark, which contains deep mutational scanning datasets that measure
protein fitness (in different contexts) for sequence variants of a wide range of proteins. This dataset contains
only single-substitution variants for the protein Anti-CRISPR protein AcrIIA4 from the organism Listeria monocytogenes. This protein has Uniprot ID: A0A247D711_LISMN. 

The DMS selection assay was described as follows: 

    Activity against spycas9 inducing an rfp reporter

It was categorised as measuring the following (general) fitness attribute: Activity. Higher scores indicate better fitness.

The source publication for this dataset is titled: 

"A deep mutational scanning platform to characterize the fitness landscape of anti-CRISPR proteins"

and can be accessed at the following DOI: 10.1101/2021.08.21.457204.

## Objective

The objective of this benchmark is to train a model that can predict the fitness of unseen single-substitution sequence variants of Anti-CRISPR protein AcrIIA4.
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
