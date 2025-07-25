
# Proteingym-DMS dataset: HSP82_YEAST_Cote-Hammarlof_2020_growth-H2O2

## Description

This dataset is part of the ProteinGym DMS benchmark, which contains deep mutational scanning datasets that measure
protein fitness (in different contexts) for sequence variants of a wide range of proteins. This dataset contains
only single-substitution variants for the protein HSP82 from the organism Saccharomyces cerevisiae. This protein has Uniprot ID: HSP82_YEAST. 

The DMS selection assay was described as follows: 

    Growth (h2o2 stress)

It was categorised as measuring the following (general) fitness attribute: OrganismalFitness. Higher scores indicate better fitness.

The source publication for this dataset is titled: 

"The Adaptive Potential of the Middle Domain of Yeast Hsp90"

and can be accessed at the following DOI: 10.1093/molbev/msaa211.

## Objective

The objective of this benchmark is to train a model that can predict the fitness of unseen single-substitution sequence variants of HSP82.
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
