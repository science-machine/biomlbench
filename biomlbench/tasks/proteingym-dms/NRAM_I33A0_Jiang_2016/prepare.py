
"""Auto-generated preparation for proteingym-dms/NRAM_I33A0_Jiang_2016."""

from pathlib import Path
import pandas as pd
import numpy as np 

def prepare(raw: Path, public: Path, private: Path) -> None:
    """
    Prepare the proteingym-dms/NRAM_I33A0_Jiang_2016 dataset from ProteinGym data.

    Args:
        raw: Directory with the DMS dataset as a CSV file
        public: Directory for public data (train.csv, test_features.csv)  
        private: Directory for private data (answers.csv)
    """

    # Load the ProteinGym DMS data files (downloaded by ProteinGymDMSDataSource)
    # The raw parameter now points to the task's raw directory containing the specific CSV
    df = pd.read_csv(raw / "NRAM_I33A0_Jiang_2016.csv")
    # Metadata is stored in the shared proteingym-dms directory
    metadata = pd.read_csv(raw.parent.parent / "DMS_substitutions.csv", index_col="DMS_id")
    fold_columns = ["fold_random_5", "fold_modulo_5", "fold_contiguous_5"]
    seq_column = "mutated_sequence"
    target_column = "DMS_score"

    wt_sequence = metadata.loc["NRAM_I33A0_Jiang_2016", "target_seq"]
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
