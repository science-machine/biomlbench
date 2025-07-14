
# Proteingym-DMS dataset: A4GRB6_PSEAI_Chen_2020

## Description

This dataset is part of the ProteinGym DMS benchmark, which contains deep mutational scanning datasets that measure
protein fitness (in different contexts) for sequence variants of a wide range of proteins. This dataset contains
only single-substitution variants for the protein Beta-lactamase VIM-2 from the organism Pseudomonas aeruginosa. This protein has Uniprot ID: A4GRB6_PSEAI. 

The DMS selection assay was described as follows: 

    Drug resistance (128/16/2.0 ug/ml ampicillin, 4.0/0.5 ug/ml cefotaxime, 0.031 ug/ml meropenem @ 25c, 37c)

It was categorised as measuring the following (general) fitness attribute: OrganismalFitness. Higher scores indicate better fitness.

The source publication for this dataset is titled: 

"Comprehensive exploration of the translocation, stability and substrate recognition requirements in VIM-2 lactamase"

and can be accessed at the following DOI: 10.7554/eLife.56707.

## Objective

The objective of this benchmark is to train a model that can predict the fitness of unseen single-substitution sequence variants of Beta-lactamase VIM-2.
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
