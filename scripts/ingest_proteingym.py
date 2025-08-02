import argparse
import tempfile
from pathlib import Path
from string import Template
from textwrap import dedent

import pandas as pd
import requests

from biomlbench.data import download_and_prepare_dataset
from biomlbench.registry import registry


def parse_args():
    parser = argparse.ArgumentParser(description="Ingest ProteinGym-DMS tasks")
    parser.add_argument(
        "--prepare",
        action="store_true",
        default=True,
        help="Also prepare datasets after creating tasks (default: True)",
    )
    parser.add_argument(
        "--no-prepare",
        dest="prepare",
        action="store_false",
        help="Skip dataset preparation (only create task structure)",
    )
    args = parser.parse_args()

    return args


def ingest_task(dms_id: str, metadata: pd.DataFrame, task_dir: Path, prepare: bool) -> str:
    """Prepares a ProteinGym DMS task."""
    script_template = Template(
        """
    \"\"\"Auto-generated preparation for proteingym-dms/$dataset_name.\"\"\"

    from pathlib import Path
    import pandas as pd
    import numpy as np 

    def prepare(raw: Path, public: Path, private: Path) -> None:
        \"\"\"
        Prepare the proteingym-dms/$dataset_name dataset from ProteinGym data.

        Args:
            raw: Directory with the DMS dataset as a CSV file
            public: Directory for public data (train.csv, test_features.csv)  
            private: Directory for private data (answers.csv)
        \"\"\"

        # Load the ProteinGym DMS data files (downloaded by ProteinGymDMSDataSource)
        # The raw parameter now points to the task's raw directory containing the specific CSV
        df = pd.read_csv(raw / "$dataset_name.csv")
        
        # Metadata is stored in the shared proteingym-dms directory
        metadata = pd.read_csv(raw.parent.parent / "DMS_substitutions.csv", index_col="DMS_id")
        fold_columns = ["fold_random_5", "fold_modulo_5", "fold_contiguous_5"]
        seq_column = "mutated_sequence"
        target_column = "DMS_score"
        
        wt_sequence = metadata.loc["$dataset_name", "target_seq"]
        wt_seq_row = pd.DataFrame({"sequence": [wt_sequence], "fitness_score": [np.nan]})
        wt_seq_row[fold_columns] = -1
        
        df["id"] = range(len(df))
        df.rename({seq_column: "sequence", target_column: "fitness_score"}, axis=1, inplace=True)
        df = df[["id", "sequence", "fitness_score"] + fold_columns]
        df.to_csv(public / "data.csv", index=False)
        
        sample_submission = pd.DataFrame({
            "id": df["id"], 
            "fitness_score_fold_random_5": [0.0] * len(df),
            "fitness_score_fold_modulo_5": [0.0] * len(df),
            "fitness_score_fold_contiguous_5": [0.0] * len(df)
        })
        sample_submission.to_csv(public / "sample_submission.csv", index=False)
        
        answers = pd.DataFrame({
            "id": df["id"],
            "fitness_score_fold_random_5": df["fitness_score"],
            "fitness_score_fold_modulo_5": df["fitness_score"],
            "fitness_score_fold_contiguous_5": df["fitness_score"]
        })
        answers.to_csv(private / "answers.csv", index=False)

        print("✅ Datasets prepared successfully!")
        print(f"Files created:")
        print(f"  Public: data.csv, sample_submission.csv")
        print(f"  Sequence column: sequence")
        print(f"  Target column: fitness_score")
        print(f"  Fold columns: {fold_columns}")
    """
    )

    config_template = Template(
        """
    id: proteingym-dms/$dataset_name
    name: $dataset_name
    task_type: fitness_prediction
    domain: protein_engineering
    difficulty: medium
    awards_medals: false
    prizes: null
    description: biomlbench/tasks/proteingym-dms/$dataset_name/description.md
    preparer: biomlbench.tasks.proteingym-dms.$dataset_name.prepare:prepare

    data_source:
      type: proteingym-dms
      benchmark_id: $dataset_name

    dataset:
      answers: proteingym-dms/$dataset_name/prepared/private/answers.csv
      sample_submission: proteingym-dms/$dataset_name/prepared/public/sample_submission.csv

    grader:
      name: proteingym-metric
      grade_fn: biomlbench.tasks.proteingym-dms.$dataset_name.grade:grade

    biomedical_metadata:
      modality: protein_sequence
      data_type: single_task
      proteingym_main_metric: spearman

    compute_requirements:
      recommended_gpu_memory_gb: 4
      estimated_runtime_minutes: 15
      max_dataset_size_gb: 2
    """
    )

    description_template = Template(
        """
    # Proteingym-DMS dataset: $dataset_name

    ## Description

    This dataset is part of the ProteinGym DMS benchmark, which contains deep mutational scanning datasets that measure
    protein fitness (in different contexts) for sequence variants of a wide range of proteins. This dataset contains
    only single-substitution variants for the protein $protein_name from the organism $organism. This protein has Uniprot ID: $uniprot_id. 

    The DMS selection assay was described as follows: 

        $selection_assay

    It was categorised as measuring the following (general) fitness attribute: $fitness_attribute. Higher scores indicate better fitness.

    The source publication for this dataset is titled: 

    "$publication_title"

    and can be accessed at the following DOI: $doi.

    ## Objective

    The objective of this benchmark is to train a model that can predict the fitness of unseen single-substitution sequence variants of $protein_name.
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
    - `sample_submission.csv`: Example submission format with ID column and fitness score column

    ## Evaluation

    Your model will be evaluated on the Spearman correlation between the predicted fitness scores and the true fitness scores for
    each of the sequences in `data.csv`.
    """
    )

    grade_template = Template(
        """
    \"\"\"Auto-generated grading for proteingym-dms/$dataset_name.\"\"\"
    import pandas as pd
    import polaris as po
    from scipy.stats import spearmanr

    def grade(submission: pd.DataFrame, answers: pd.DataFrame) -> float:
        \"\"\"Grade using spearman correlation between true and predicted fitness scores.\"\"\"

        # Extract predictions for each fold
        fold_random_5 = submission[submission.columns[1]].tolist()
        fold_modulo_5 = submission[submission.columns[2]].tolist()
        fold_contiguous_5 = submission[submission.columns[3]].tolist()
        answers_fold_random_5 = answers[answers.columns[1]].tolist()
        answers_fold_modulo_5 = answers[answers.columns[2]].tolist()
        answers_fold_contiguous_5 = answers[answers.columns[3]].tolist()

        # Calculate Spearman correlation for each split
        correlation_fold_random_5, _ = spearmanr(fold_random_5, answers_fold_random_5)
        correlation_fold_modulo_5, _ = spearmanr(fold_modulo_5, answers_fold_modulo_5)
        correlation_fold_contiguous_5, _ = spearmanr(fold_contiguous_5, answers_fold_contiguous_5)

        # Return the average correlation across all splits
        return (correlation_fold_random_5 + correlation_fold_modulo_5 + correlation_fold_contiguous_5) / 3
    """
    )

    (task_dir / dms_id).mkdir(parents=True, exist_ok=True)
    with open(task_dir / dms_id / "prepare.py", "w") as f:
        f.write(dedent(script_template.substitute(dataset_name=dms_id)))
    with open(task_dir / dms_id / "config.yaml", "w") as f:
        f.write(dedent(config_template.substitute(dataset_name=dms_id)))
    with open(task_dir / dms_id / "description.md", "w") as f:
        f.write(
            dedent(
                description_template.substitute(
                    dataset_name=dms_id,
                    protein_name=metadata.loc[dms_id, "molecule_name"],
                    organism=metadata.loc[dms_id, "source_organism"],
                    uniprot_id=metadata.loc[dms_id, "UniProt_ID"],
                    selection_assay=(
                        metadata.loc[dms_id, "selection_assay"].capitalize()
                        if isinstance(metadata.loc[dms_id, "selection_assay"], str)
                        else "(Not provided)"
                    ),
                    fitness_attribute=metadata.loc[dms_id, "coarse_selection_type"],
                    publication_title=metadata.loc[dms_id, "title"],
                    doi=metadata.loc[dms_id, "jo"],
                )
            )
        )
    with open(task_dir / dms_id / "grade.py", "w") as f:
        f.write(dedent(grade_template.substitute(dataset_name=dms_id)))

    if prepare:
        print(f"🔄 Preparing data for proteingym/{dms_id}...")
        task = registry.get_task(f"proteingym-dms/{dms_id}")
        download_and_prepare_dataset(task)

    return dms_id


def main():
    """
    Ingest the proteingym-dms dataset into the biomlbench format.

    This script will create a directory for each DMS dataset in ProteinGym which contains only single-substitution variants.
    """

    TASK_DIR = Path("biomlbench/tasks/proteingym-dms")

    # URL for the metadata CSV
    DMS_SUBSTITUTIONS_URL = "https://zenodo.org/records/15293562/files/DMS_substitutions.csv"

    args = parse_args()

    # Download the metadata CSV to a temporary file and load
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmpfile:
        response = requests.get(DMS_SUBSTITUTIONS_URL)
        response.raise_for_status()
        tmpfile.write(response.content)
        tmp_csv_path = tmpfile.name

    dms_substitutions = pd.read_csv(tmp_csv_path).set_index("DMS_id")
    dms_single_substitutions = dms_substitutions[
        dms_substitutions["includes_multiple_mutants"] == False
    ].copy()
    dms_ids = dms_single_substitutions.index.tolist()

    for dms_id in dms_ids:
        try:
            dms_id = ingest_task(dms_id, dms_single_substitutions, TASK_DIR, args.prepare)
            status = "✅ Created and prepared" if args.prepare else "✅ Created"
            print(f"{status} proteingym-dms/{dms_id.replace('/', '-').lower()}")
        except Exception as e:
            print(f"❌ Failed to create {dms_id}: {e}")


if __name__ == "__main__":
    main()
