
# Proteingym-DMS dataset: HEM3_HUMAN_Loggerenberg_2023

## Description

This dataset is part of the ProteinGym DMS benchmark, which contains deep mutational scanning datasets that measure
protein fitness (in different contexts) for sequence variants of a wide range of proteins. This dataset contains
only single-substitution variants for the protein hydroxymethylbilane synthase from the organism Homo sapiens. This protein has Uniprot ID: HEM3_HUMAN. 

The DMS selection assay was described as follows: 

    Activity

It was categorised as measuring the following (general) fitness attribute: Activity. Higher scores indicate better fitness.

The source publication for this dataset is titled: 

"Systematically testing human HMBS missense variants to reveal mechanism and pathogenic variation"

and can be accessed at the following DOI: 10.1101/2023.02.06.527353.

## Objective

The objective of this benchmark is to train a model that can predict the fitness of unseen single-substitution sequence variants of hydroxymethylbilane synthase.
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
