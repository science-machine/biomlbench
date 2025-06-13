"""
Caco-2 Cell Permeability Prediction Task

A drug discovery task for predicting intestinal permeability using Caco-2 cell assays.
Based on the TDC/Polaris benchmark.

Task Details:
- Input: SMILES molecular representations
- Output: Caco-2 permeability values (regression)
- Metric: Mean Absolute Error (MAE)
- Domain: Drug Discovery/ADMET
- Clinical Relevance: Oral drug absorption prediction

The task requires predicting how well drug molecules will permeate
through intestinal tissue, which is critical for oral bioavailability.
"""

from .grade import grade, evaluate_with_polaris
from .prepare import prepare

__all__ = ["grade", "prepare", "evaluate_with_polaris"] 