"""
Task-specific baselines for histopathologic-cancer-detection medical imaging task.

These baselines are specific to medical imaging tasks and demonstrate
different approaches to image classification.
"""

import pandas as pd
import numpy as np

from biomlbench.baselines import BaselineAgent, BaselineFactory
from biomlbench.registry import Task
from biomlbench.utils import get_logger

logger = get_logger(__name__)


class ImageClassificationBaseline(BaselineAgent):
    """Simple image classification baseline using basic CNN."""
    
    def __init__(self, seed: int = 42):
        super().__init__("cnn", seed)
    
    def generate_submission(self, task: Task) -> pd.DataFrame:
        """Placeholder for CNN-based image classification."""
        # Note: This is a conceptual example - actual implementation would
        # require loading images, training CNN, etc.
        
        # Load sample submission to get format
        sample_submission = pd.read_csv(task.sample_submission)
        
        # For demonstration: predict based on random probabilities
        # In real implementation, this would be CNN predictions
        n_predictions = len(sample_submission)
        predictions = np.random.random(n_predictions)  # Random probabilities
        
        # Create submission
        submission = sample_submission.copy()
        prediction_col = [col for col in submission.columns if col.lower() != 'id'][0]
        submission[prediction_col] = predictions
        
        logger.info(f"Generated {n_predictions} CNN baseline predictions (placeholder)")
        return submission


class ResNetBaseline(BaselineAgent):
    """ResNet-based baseline for medical imaging."""
    
    def __init__(self, seed: int = 42):
        super().__init__("resnet", seed)
    
    def generate_submission(self, task: Task) -> pd.DataFrame:
        """Placeholder for ResNet-based classification."""
        # Load sample submission
        sample_submission = pd.read_csv(task.sample_submission)
        
        # Placeholder predictions
        n_predictions = len(sample_submission)
        predictions = np.random.random(n_predictions) * 0.8 + 0.1  # Slightly different distribution
        
        # Create submission
        submission = sample_submission.copy()
        prediction_col = [col for col in submission.columns if col.lower() != 'id'][0]
        submission[prediction_col] = predictions
        
        logger.info(f"Generated {n_predictions} ResNet baseline predictions (placeholder)")
        return submission


def register_baselines():
    """Register histopathologic-cancer-detection specific baselines with the factory."""
    BaselineFactory.register_baseline("histopathologic-cancer-detection", "cnn", ImageClassificationBaseline)
    BaselineFactory.register_baseline("histopathologic-cancer-detection", "resnet", ResNetBaseline)
    logger.info("Registered medical imaging baselines for histopathologic-cancer-detection task") 