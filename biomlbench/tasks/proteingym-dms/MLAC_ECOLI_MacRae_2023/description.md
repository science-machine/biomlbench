
# Proteingym-DMS dataset: MLAC_ECOLI_MacRae_2023

## Description

This dataset is part of the ProteinGym DMS benchmark, which contains deep mutational scanning datasets that measure
protein fitness (in different contexts) for sequence variants of a wide range of proteins. This dataset contains
only single-substitution variants for the protein MlaC lipid transporter from the organism Escherichia coli. This protein has Uniprot ID: MLAC_ECOLI. 

The DMS selection assay was described as follows: 

    Cell growth in ∆mlac and selective medium

It was categorised as measuring the following (general) fitness attribute: OrganismalFitness. Higher scores indicate better fitness.

The source publication for this dataset is titled: 

"Protein-protein interactions in the Mla lipid transport system probed by computational structure prediction and deep mutational scanning"

and can be accessed at the following DOI: 10.1016/j.jbc.2023.104744.

## Objective

The objective of this benchmark is to train a model that can predict the fitness of unseen single-substitution sequence variants of MlaC lipid transporter.
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
