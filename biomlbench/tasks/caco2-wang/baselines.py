"""
Task-specific baselines for caco2-wang molecular permeability prediction.

These baselines are specific to molecular tasks with SMILES strings and
demonstrate different approaches to molecular property prediction.
"""

from typing import Optional
import pandas as pd
import numpy as np

from biomlbench.baselines import BaselineAgent, BaselineFactory
from biomlbench.registry import Task
from biomlbench.utils import get_logger

logger = get_logger(__name__)


class MolecularLinearBaseline(BaselineAgent):
    """Linear regression baseline using simple molecular features."""
    
    def __init__(self, seed: int = 42):
        super().__init__("linear", seed)
    
    def generate_submission(self, task: Task) -> pd.DataFrame:
        """Use linear regression with simple molecular features."""
        from sklearn.linear_model import LinearRegression
        from sklearn.preprocessing import StandardScaler
        from sklearn.impute import SimpleImputer
        
        # Load data
        train_df = pd.read_csv(task.public_dir / "train.csv")
        test_df = pd.read_csv(task.public_dir / "test_features.csv")
        sample_submission = pd.read_csv(task.sample_submission)
        
        # Extract features and target
        target_cols = [col for col in train_df.columns if col not in ['id', 'smiles']]
        target_col = target_cols[0]
        
        # Simple molecular features (string-based, no chemistry libraries needed)
        def extract_simple_features(smiles_list):
            features = []
            for smiles in smiles_list:
                feat = {
                    'length': len(smiles),
                    'num_C': smiles.count('C'),
                    'num_N': smiles.count('N'),
                    'num_O': smiles.count('O'),
                    'num_S': smiles.count('S'),
                    'num_rings': smiles.count('c'),  # Approximate
                    'num_branches': smiles.count('('),
                    'num_double_bonds': smiles.count('='),
                    'num_triple_bonds': smiles.count('#'),
                    'complexity': len(set(smiles)),
                }
                features.append(feat)
            return pd.DataFrame(features)
        
        # Extract features
        X_train = extract_simple_features(train_df['smiles'])
        X_test = extract_simple_features(test_df['smiles'])
        y_train = train_df[target_col]
        
        # Handle missing values and scale
        imputer = SimpleImputer(strategy='median')
        scaler = StandardScaler()
        
        X_train_processed = scaler.fit_transform(imputer.fit_transform(X_train))
        X_test_processed = scaler.transform(imputer.transform(X_test))
        
        # Train model
        model = LinearRegression()
        model.fit(X_train_processed, y_train)
        
        # Make predictions
        predictions = model.predict(X_test_processed)
        
        # Create submission
        submission = sample_submission.copy()
        prediction_col = [col for col in submission.columns if col != 'id'][0]
        submission[prediction_col] = predictions
        
        logger.info(f"Generated {len(predictions)} molecular linear regression predictions")
        return submission


class MolecularRandomForestBaseline(BaselineAgent):
    """Random Forest baseline with enhanced molecular features."""
    
    def __init__(self, seed: int = 42, n_estimators: int = 100, max_depth: Optional[int] = None):
        super().__init__("rf", seed)
        self.n_estimators = n_estimators
        self.max_depth = max_depth
    
    def generate_submission(self, task: Task) -> pd.DataFrame:
        """Use Random Forest with enhanced molecular features."""
        from sklearn.ensemble import RandomForestRegressor
        from sklearn.impute import SimpleImputer
        
        # Load data
        train_df = pd.read_csv(task.public_dir / "train.csv")
        test_df = pd.read_csv(task.public_dir / "test_features.csv")
        sample_submission = pd.read_csv(task.sample_submission)
        
        # Extract features and target
        target_cols = [col for col in train_df.columns if col not in ['id', 'smiles']]
        target_col = target_cols[0]
        
        # Enhanced molecular features
        def extract_enhanced_features(smiles_list):
            features = []
            for smiles in smiles_list:
                feat = {
                    'length': len(smiles),
                    'num_C': smiles.count('C'),
                    'num_N': smiles.count('N'),
                    'num_O': smiles.count('O'),
                    'num_S': smiles.count('S'),
                    'num_P': smiles.count('P'),
                    'num_F': smiles.count('F'),
                    'num_Cl': smiles.count('Cl'),
                    'num_Br': smiles.count('Br'),
                    'num_I': smiles.count('I'),
                    'num_rings': smiles.count('c'),
                    'num_branches': smiles.count('('),
                    'num_double_bonds': smiles.count('='),
                    'num_triple_bonds': smiles.count('#'),
                    'complexity': len(set(smiles)),
                    'aromatic_ratio': smiles.count('c') / max(len(smiles), 1),
                    'heteroatom_ratio': (smiles.count('N') + smiles.count('O') + smiles.count('S')) / max(len(smiles), 1),
                    'branch_ratio': smiles.count('(') / max(len(smiles), 1),
                }
                features.append(feat)
            return pd.DataFrame(features)
        
        # Extract features
        X_train = extract_enhanced_features(train_df['smiles'])
        X_test = extract_enhanced_features(test_df['smiles'])
        y_train = train_df[target_col]
        
        # Handle missing values
        imputer = SimpleImputer(strategy='median')
        X_train_processed = imputer.fit_transform(X_train)
        X_test_processed = imputer.transform(X_test)
        
        # Train model
        model = RandomForestRegressor(
            n_estimators=self.n_estimators,
            max_depth=self.max_depth,
            random_state=self.seed,
            n_jobs=-1
        )
        model.fit(X_train_processed, y_train)
        
        # Make predictions
        predictions = model.predict(X_test_processed)
        
        # Create submission
        submission = sample_submission.copy()
        prediction_col = [col for col in submission.columns if col != 'id'][0]
        submission[prediction_col] = predictions
        
        logger.info(f"Generated {len(predictions)} molecular Random Forest predictions (n_estimators={self.n_estimators})")
        return submission


class MolecularFingerprintBaseline(BaselineAgent):
    """Baseline using molecular fingerprint features (if RDKit available)."""
    
    def __init__(self, seed: int = 42):
        super().__init__("fingerprint", seed)
    
    def generate_submission(self, task: Task) -> pd.DataFrame:
        """Use molecular fingerprints with Random Forest."""
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors
            from sklearn.ensemble import RandomForestRegressor
            from sklearn.impute import SimpleImputer
        except ImportError:
            raise ImportError("RDKit is required for fingerprint baseline. Install with: conda install -c conda-forge rdkit")
        
        # Load data
        train_df = pd.read_csv(task.public_dir / "train.csv")
        test_df = pd.read_csv(task.public_dir / "test_features.csv")
        sample_submission = pd.read_csv(task.sample_submission)
        
        # Extract features and target
        target_cols = [col for col in train_df.columns if col not in ['id', 'smiles']]
        target_col = target_cols[0]
        
        # RDKit molecular descriptors
        def extract_rdkit_features(smiles_list):
            features = []
            for smiles in smiles_list:
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is not None:
                        feat = {
                            'mw': Descriptors.MolWt(mol),
                            'logp': Descriptors.MolLogP(mol),
                            'hbd': Descriptors.NumHDonors(mol),
                            'hba': Descriptors.NumHAcceptors(mol),
                            'tpsa': Descriptors.TPSA(mol),
                            'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
                            'aromatic_rings': Descriptors.NumAromaticRings(mol),
                            'heavy_atoms': Descriptors.HeavyAtomCount(mol),
                        }
                    else:
                        # Fallback for invalid SMILES
                        feat = {k: 0.0 for k in ['mw', 'logp', 'hbd', 'hba', 'tpsa', 'rotatable_bonds', 'aromatic_rings', 'heavy_atoms']}
                except:
                    # Fallback for any errors
                    feat = {k: 0.0 for k in ['mw', 'logp', 'hbd', 'hba', 'tpsa', 'rotatable_bonds', 'aromatic_rings', 'heavy_atoms']}
                features.append(feat)
            return pd.DataFrame(features)
        
        # Extract features
        X_train = extract_rdkit_features(train_df['smiles'])
        X_test = extract_rdkit_features(test_df['smiles'])
        y_train = train_df[target_col]
        
        # Handle missing values
        imputer = SimpleImputer(strategy='median')
        X_train_processed = imputer.fit_transform(X_train)
        X_test_processed = imputer.transform(X_test)
        
        # Train model
        model = RandomForestRegressor(
            n_estimators=100,
            max_depth=10,
            random_state=self.seed,
            n_jobs=-1
        )
        model.fit(X_train_processed, y_train)
        
        # Make predictions
        predictions = model.predict(X_test_processed)
        
        # Create submission
        submission = sample_submission.copy()
        prediction_col = [col for col in submission.columns if col != 'id'][0]
        submission[prediction_col] = predictions
        
        logger.info(f"Generated {len(predictions)} molecular fingerprint predictions")
        return submission


def register_baselines():
    """Register caco2-wang specific baselines with the factory."""
    BaselineFactory.register_baseline("caco2-wang", "linear", MolecularLinearBaseline)
    BaselineFactory.register_baseline("caco2-wang", "rf", MolecularRandomForestBaseline)
    BaselineFactory.register_baseline("caco2-wang", "fingerprint", MolecularFingerprintBaseline)
    logger.info("Registered molecular baselines for caco2-wang task") 