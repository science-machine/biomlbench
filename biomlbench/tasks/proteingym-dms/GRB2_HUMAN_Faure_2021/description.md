
# Proteingym-DMS dataset: GRB2_HUMAN_Faure_2021

## Description

This dataset is part of the ProteinGym DMS benchmark, which contains deep mutational scanning datasets that measure
protein fitness (in different contexts) for sequence variants of a wide range of proteins. This dataset contains
multi-substitution variants for the protein GRB2-SH3 from the organism Homo sapiens. This protein has Uniprot ID: GRB2_HUMAN. 

The DMS selection assay was described as follows: 

    Yeast growth

It was categorised as measuring the following (general) fitness attribute: OrganismalFitness. Higher scores indicate better fitness.

The source publication for this dataset is titled: 

"Mapping the energetic and allosteric landscapes of protein binding domains"

and can be accessed at the following DOI: 10.1038/s41586-022-04586-4.

## Objective

The objective of this benchmark is to train a model that can predict the fitness of unseen single-substitution sequence variants of GRB2-SH3.
To train your model, you will use 5-fold cross-validation on the sequences and fitness scores defined in the `data.csv` file. 

You will use the `fold_random_5`, `fold_modulo_5`, and `fold_contiguous_5` columns to split the data into training and test sets.
Each of these columns contains integer values from 0 to 4, which indicate the fold of the sequence in the corresponding 5-fold cross-validation split.

When predicting the fitness score for a given sequence, **you must use a model trained only on sequences from other folds**.
For example, to predict the fitness score for sequences in fold 0 for the `fold_random_5` split, you must use a model trained
only on the sequences with `fold_random_5` values other than 0.

You must repeat this process for each of the five folds in `fold_random_5` (so that all sequences in `data.csv` 
receive a predicted score). Then repeat the process separately for the other two cross-validation split columns
`fold_modulo_5` and `fold_contiguous_5`. Hence, each sequence should have three predicted fitness scores,
corresponding to the prediction for that sequence under models trained on the three different cross-validation splits.

Overall, your training and inference pseudocode loop should look like this:

```python
import pandas as pd
data = pd.read_csv("data.csv")
wt_sequence = data.iloc[0]["sequence"]

# remove the wild-type sequence from the data
data = data.iloc[1:]

## define your model here ##
model = ...

# initialize a dataframe to store the predictions
predictions = pd.DataFrame(columns=["id", "fold_random_5", "fold_modulo_5", "fold_contiguous_5"], index=data.index)
predictions["id"] = data["id"]

# loop over the different cross-validation splits
fold_columns = ["fold_random_5", "fold_modulo_5", "fold_contiguous_5"]
for column in fold_columns:  # different splits
    for fold in range(5):  # different folds per-split
        fold_mask = data[column] == fold
        train_data = data[~fold_mask]  # train on all folds except the current one
        test_data = data[fold_mask]  # test on the current fold

        # train the model **from scratch** on the training set 
        trained_model = model.fit(train_data["sequence"], train_data["fitness_score"]) 

        # predict the fitness score for the sequences in the current fold 
        predictions.loc[fold_mask, f"fitness_score_{column}"] = trained_model.predict(test_data["sequence"])
```

Hence, the output data frame should contain four columns:
- `id`: The ID of the sequence 
- `fitness_score_fold_random_5`: The predicted fitness score for that sequence predicted by the model trained on the 
  cross-validation split defined by the `fold_random_5` column
- `fitness_score_fold_modulo_5`: The predicted fitness score for that sequence predicted by the model trained on the 
  cross-validation split defined by the `fold_modulo_5` column
- `fitness_score_fold_contiguous_5`: The predicted fitness score for that sequence predicted by the model trained on the 
  cross-validation split defined by the `fold_contiguous_5` column

## Data Format

The `data.csv` file contains the following columns:
- `id`: The index of the sequence
- `sequence`: The amino acid sequence of the variant
- `fitness_score`: The fitness score of the variant
- `fold_random_5`: The fold of the variant (0-4) in the "random" 5-fold cross-validation split
- `fold_modulo_5`: The fold of the variant (0-4) in the "modulo" 5-fold cross-validation split
- `fold_contiguous_5`: The fold of the variant (0-4) in the "contiguous" 5-fold cross-validation split

The first row of the CSV file contains the wild-type sequence in the `sequence` field and missing values for the other columns.

## Files

- `data.csv`: File with sequences and fitness scores
- `sample_submission.csv`: Example submission format with ID column and the three fitness score columns

## Evaluation

Your model will be evaluated on the Spearman correlation between the predicted fitness scores and the true fitness scores for
each of the sequences in `data.csv`.
