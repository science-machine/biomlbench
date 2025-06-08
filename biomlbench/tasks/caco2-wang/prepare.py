"""
Data preparation for Caco-2 Wang permeability prediction task.

This script processes Polaris data downloaded by the data source system
and converts it to BioML-bench format.
"""

from pathlib import Path

import pandas as pd


def prepare(raw: Path, public: Path, private: Path) -> None:
    """
    Prepare the Caco-2 permeability prediction dataset from Polaris data.
    
    Args:
        raw: Directory with raw Polaris data (downloaded by PolarisDataSource)
        public: Directory for public data (train.csv, test_features.csv)  
        private: Directory for private data (answers.csv)
    """
    
    # Load the Polaris data files (downloaded by PolarisDataSource)
    train_df = pd.read_csv(raw / 'polaris_train_data.csv')
    test_df = pd.read_csv(raw / 'polaris_test_data.csv')
    
    print(f"Loaded Polaris data: {len(train_df)} train, {len(test_df)} test samples")
    
    # Create public training data (train.csv)
    # Format: SMILES, target
    train_public = train_df.copy()
    train_public = train_public.rename(columns={'Drug': 'smiles', 'Y': 'caco2_permeability'})
    train_public.to_csv(public / 'train.csv', index=False)
    
    # Create test features (test_features.csv) 
    # Format: id, SMILES (no targets)
    test_features = test_df.copy()
    test_features = test_features.rename(columns={'Drug': 'smiles'})
    test_features['id'] = range(len(test_features))  # Add ID column
    test_features = test_features[['id', 'smiles']]  # Reorder columns
    test_features.to_csv(public / 'test_features.csv', index=False)
    
    # Create sample submission (sample_submission.csv)
    sample_submission = pd.DataFrame({
        'id': range(len(test_features)),
        'caco2_permeability': [0.0] * len(test_features)  # Dummy predictions
    })
    sample_submission.to_csv(public / 'sample_submission.csv', index=False)
    
    # Create private answers file for evaluation
    answers = pd.DataFrame({
        'id': range(len(test_features)),
        'caco2_permeability': test_df['Y'].values  # Use actual test labels
    })
    answers.to_csv(private / 'answers.csv', index=False)
    
    print(f"Dataset prepared successfully!")
    print(f"Training examples: {len(train_public)}")
    print(f"Test examples: {len(test_features)}")
    print(f"Files created:")
    print(f"  Public: train.csv, test_features.csv, sample_submission.csv")
    print(f"  Private: answers.csv") 