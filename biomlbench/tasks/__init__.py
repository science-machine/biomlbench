"""
BioML-bench tasks module.

This module contains all biomedical task implementations for the BioML-bench framework.
Tasks are organized by domain and include:

- Medical Imaging: histopathologic-cancer-detection, medical-segmentation, etc.
- Protein Engineering: proteingym-fitness, structure-prediction, etc.  
- Drug Discovery: molecular-property-prediction, admet-prediction, etc.
- Clinical Biomarkers: biomarker-discovery, clinical-prediction, etc.
"""

# Utilities for task implementations
from . import utils

__all__ = ["utils"] 