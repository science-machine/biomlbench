
# Proteingym-DMS dataset: SPIKE_SARS2_Starr_2020_binding

## Description

This dataset is part of the ProteinGym DMS benchmark, which contains deep mutational scanning datasets that measure
protein fitness (in different contexts) for sequence variants of a wide range of proteins. This dataset contains
only single-substitution variants for the protein Spike RBD from the organism Severe acute respiratory syndrome coronavirus 2 (2019-nCoV) (SARS-CoV-2). This protein has Uniprot ID: SPIKE_SARS2. 

The DMS selection assay was described as follows: 

    Ace2 binding

It was categorised as measuring the following (general) fitness attribute: Binding. Higher scores indicate better fitness.

The source publication for this dataset is titled: 

"Deep Mutational Scanning of SARS-CoV-2 Receptor Binding Domain Reveals Constraints on Folding and ACE2 Binding"

and can be accessed at the following DOI: 10.1016/j.cell.2020.08.012.

## Objective

The objective of this benchmark is to train a model that can predict the fitness of unseen single-substitution sequence variants of Spike RBD.
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
