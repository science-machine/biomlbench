import argparse
import tempfile
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import TimeoutError as FuturesTimeoutError
from concurrent.futures import as_completed
from pathlib import Path
from string import Template
from textwrap import dedent

import pandas as pd
import requests

from biomlbench.data import download_and_prepare_datasets
from biomlbench.registry import registry


def parse_args():
    parser = argparse.ArgumentParser(description="Ingest Polaris tasks")
    parser.add_argument(
        "--workers", type=int, default=4, help="Number of workers to use for parallel processing"
    )
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


def ingest_task(dms_id: str, metadata: pd.DataFrame, task_dir: Path, prepare: bool):
    """Prepares a ProteinGym DMS task."""
    script_template = Template(
        """
    \"\"\"Auto-generated preparation for proteingym-dms/$dataset_name.\"\"\"

    from pathlib import Path
    import pandas as pd

    def prepare(raw: Path, prepared: Path) -> None:
        \"\"\"
        Prepare the proteingym-dms/$dataset_name dataset from ProteinGym data.

        Args:
            raw: Directory with the DMS dataset as a CSV file
            prepared: Directory for prepared data (train.csv, test_features.csv, sample_submission.csv, answers.csv)
        \"\"\"

        # Load the ProteinGym DMS data files (downloaded by ProteinGymDMSDataSource)
        df = pd.read_csv(raw / "cv_folds_singles_substitutions" / "$dataset_name.csv")
        fold_columns = ["fold_random_5", "fold_modulo_5", "fold_contiguous_5"]
        seq_column = "mutated_sequence"
        target_column = "DMS_score"

        for fold_column in fold_columns:
            fold_name = fold_column.rstrip('_5')
            for fold, fold_df in df.groupby(fold_column):
                (prepared / f"{fold_name}_{fold}").mkdir(parents=True, exist_ok=True)
                public = prepared / f"{fold_name}_{fold}" / "public"
                private = prepared / f"{fold_name}_{fold}" / "private"
                public.mkdir(parents=True, exist_ok=True)
                private.mkdir(parents=True, exist_ok=True)

                train_df = df[df[fold_column] != fold].copy()
                train_df["id"] = range(len(train_df))
                train_df.rename({seq_column: "sequence", target_column: "fitness_score"}, axis=1, inplace=True)

                fold_df["id"] = range(len(fold_df))
                fold_df.rename({seq_column: "sequence", target_column: "fitness_score"}, axis=1, inplace=True)

                train_df[["id", "sequence", "fitness_score"]].to_csv(public / "train.csv", index=False)
                fold_df[["id", "sequence"]].to_csv(public / "test_features.csv", index=False)
                fold_df[["id", "fitness_score"]].to_csv(private / "answers.csv", index=False)

                sample_submission = pd.DataFrame({
                    "id": fold_df["id"],
                    "fitness_score": [0.0] * len(fold_df)
                })
                sample_submission.to_csv(public / "sample_submission.csv", index=False)

        print("✅ Datasets prepared successfully!")
        print(f"Files created:")
        print(f"  Public: train.csv, test_features.csv, sample_submission.csv")
        print(f"  Private: answers.csv")
        print(f"  Sequence column: sequence")
        print(f"  Target column: fitness_score")
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
      type: proteingym
      benchmark_id: $dataset_name

    datasets:
      fold_random_0:
        answers: answers.csv
        sample_submission: sample_submission.csv
      fold_random_1:
        answers: answers.csv
        sample_submission: sample_submission.csv
      fold_random_2:
        answers: answers.csv
        sample_submission: sample_submission.csv
      fold_random_3:
        answers: answers.csv
        sample_submission: sample_submission.csv
      fold_random_4:
        answers: answers.csv
        sample_submission: sample_submission.csv
      fold_modulo_0:
        answers: answers.csv
        sample_submission: sample_submission.csv
      fold_modulo_1:
        answers: answers.csv
        sample_submission: sample_submission.csv
      fold_modulo_2:
        answers: answers.csv
        sample_submission: sample_submission.csv
      fold_modulo_3:
        answers: answers.csv
        sample_submission: sample_submission.csv
      fold_modulo_4:
        answers: answers.csv
        sample_submission: sample_submission.csv
      fold_contiguous_0:
        answers: answers.csv
        sample_submission: sample_submission.csv
      fold_contiguous_1:
        answers: answers.csv
        sample_submission: sample_submission.csv
      fold_contiguous_2:
        answers: answers.csv
        sample_submission: sample_submission.csv
      fold_contiguous_3:
        answers: answers.csv
        sample_submission: sample_submission.csv
      fold_contiguous_4:
        answers: answers.csv
        sample_submission: sample_submission.csv

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

        # Extract predictions (second column)
        predictions = submission[submission.columns[1]].tolist()
        answers = answers[answers.columns[1]].tolist()

        # Calculate Spearman correlation
        correlation, _ = spearmanr(predictions, answers)

        return correlation
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
        download_and_prepare_datasets(task)


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
            ingest_task(dms_id, dms_single_substitutions, TASK_DIR, args.prepare)
            status = "✅ Created and prepared" if args.prepare else "✅ Created"
            print(f"{status} proteingym-dms/{dms_id.replace('/', '-').lower()}")
        except Exception as e:
            print(f"❌ Failed to create {dms_id}: {e}")


if __name__ == "__main__":
    main()
