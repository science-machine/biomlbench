
"""Auto-generated preparation for proteingym-dms/NCAP_I34A1_Doud_2015."""

from pathlib import Path
import pandas as pd
import numpy as np 

def prepare(raw: Path, prepared: Path) -> None:
    """
    Prepare the proteingym-dms/NCAP_I34A1_Doud_2015 dataset from ProteinGym data.

    Args:
        raw: Directory with the DMS dataset as a CSV file
        prepared: Directory for prepared data (train.csv, test_features.csv, sample_submission.csv, answers.csv)
    """

    # Load the ProteinGym DMS data files (downloaded by ProteinGymDMSDataSource)
    df = pd.read_csv(raw / "cv_folds_singles_substitutions" / "NCAP_I34A1_Doud_2015.csv")
    metadata = pd.read_csv(raw / "DMS_substitutions.csv", index_col="DMS_id")
    wt_sequence = metadata.loc["NCAP_I34A1_Doud_2015", "target_seq"]
    wt_seq_row = pd.DataFrame({"sequence": [wt_sequence], "fitness_score": [np.nan]})
    fold_columns = ["fold_random_5", "fold_modulo_5", "fold_contiguous_5"]
    seq_column = "mutated_sequence"
    target_column = "DMS_score"

    for fold_column in fold_columns:
        fold_name = fold_column.rstrip('_5')
        for fold, test_df in df.groupby(fold_column):
            (prepared / f"{fold_name}_{fold}").mkdir(parents=True, exist_ok=True)
            public = prepared / f"{fold_name}_{fold}" / "public"
            private = prepared / f"{fold_name}_{fold}" / "private"
            public.mkdir(parents=True, exist_ok=True)
            private.mkdir(parents=True, exist_ok=True)

            train_df = df[df[fold_column] != fold].copy()
            train_df.rename({seq_column: "sequence", target_column: "fitness_score"}, axis=1, inplace=True)
            train_df = train_df[["sequence", "fitness_score"]]
            train_df = pd.concat([wt_seq_row, train_df], axis=0, ignore_index=True)
            train_df["id"] = range(len(train_df))

            test_df["id"] = range(len(test_df))
            test_df.rename({seq_column: "sequence", target_column: "fitness_score"}, axis=1, inplace=True)

            train_df[["id", "sequence", "fitness_score"]].to_csv(public / "train.csv", index=False)
            test_df[["id", "sequence"]].to_csv(public / "test_features.csv", index=False)
            test_df[["id", "fitness_score"]].to_csv(private / "answers.csv", index=False)

            sample_submission = pd.DataFrame({
                "id": test_df["id"],
                "fitness_score": [0.0] * len(test_df)
            })
            sample_submission.to_csv(public / "sample_submission.csv", index=False)

    print("✅ Datasets prepared successfully!")
    print(f"Files created:")
    print(f"  Public: train.csv, test_features.csv, sample_submission.csv")
    print(f"  Private: answers.csv")
    print(f"  Sequence column: sequence")
    print(f"  Target column: fitness_score")
