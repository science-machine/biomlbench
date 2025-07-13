"""Auto-generated preparation for biogen/adme-fang-hppb-reg-v1."""
from pathlib import Path
import pandas as pd

def prepare(raw: Path, prepared: Path) -> None:
    """
    Prepare biogen-adme-fang-hppb-reg-v1 dataset from Polaris data downloaded by PolarisDataSource.
    
    Args:
        raw: Directory with raw data
        prepared: Directory for prepared data - will create 0/public and 0/private subdirectories
    """
    
    # Create fold 0 directories (for single fold, we use 0 as the default)
    dataset_dir = prepared / "0"
    public_dir = dataset_dir / "public"
    private_dir = dataset_dir / "private"
    
    public_dir.mkdir(parents=True, exist_ok=True)
    private_dir.mkdir(parents=True, exist_ok=True)
    
    # Load the Polaris data files (downloaded by PolarisDataSource)
    train_df = pd.read_csv(raw / 'polaris_train_data.csv')
    test_df = pd.read_csv(raw / 'polaris_test_data.csv')
    
    print(f"Loaded Polaris data: {len(train_df)} train, {len(test_df)} test samples")
    print(f"Columns: {list(train_df.columns)}")
    
    # First column is ALWAYS the feature/molecule column
    molecule_col = train_df.columns[0]
    
    # Load benchmark to get target metadata
    import polaris as po
    benchmark = po.load_benchmark("biogen/adme-fang-hppb-reg-v1")
    polaris_target_cols = list(benchmark.target_cols)
    
    # Use first target if multiple targets exist
    if len(polaris_target_cols) > 1:
        target_col = polaris_target_cols[0]
        print(f"Multiple targets found: {polaris_target_cols}")
        print(f"Using first target only: {target_col}")
    elif len(polaris_target_cols) == 1:
        target_col = polaris_target_cols[0]
    else:
        raise ValueError(f"No targets found in benchmark biogen/adme-fang-hppb-reg-v1")
    
    print(f"Using molecule column: {molecule_col}")
    print(f"Using target column: {target_col}")
    
    # Create public training data (train.csv)
    train_public = train_df[[molecule_col, target_col]].copy()
    train_public.to_csv(public_dir / 'train.csv', index=False)
    
    # Create test features (test_features.csv) - no targets
    test_features = test_df[[molecule_col]].copy()
    test_features['id'] = range(len(test_features))
    test_features = test_features[['id', molecule_col]]
    test_features.to_csv(public_dir / 'test_features.csv', index=False)
    
    # Create sample submission (sample_submission.csv)
    sample_submission = pd.DataFrame({
        'id': range(len(test_features)),
        target_col: [0.0] * len(test_features)
    })
    sample_submission.to_csv(public_dir / 'sample_submission.csv', index=False)
    
    # Create private answers file for evaluation
    answers = pd.DataFrame({
        'id': range(len(test_features)),
        target_col: test_df[target_col].values
    })
    answers.to_csv(private_dir / 'answers.csv', index=False)
    
    print("✅ Dataset prepared successfully!")
    print(f"Molecule column: {molecule_col}")
    print(f"Target column: {target_col}")
    print(f"Files created:")
    print(f"  Public: train.csv, test_features.csv, sample_submission.csv")
    print(f"  Private: answers.csv")
