"""
Histopathologic Cancer Detection Task

A medical imaging task for detecting metastatic cancer in histopathology patches.
Based on the PatchCamelyon (PCam) dataset.

Task Details:
- Input: 96x96 pixel histopathology images  
- Output: Binary classification (cancer/no cancer)
- Metric: Area Under ROC Curve (AUC-ROC)
- Domain: Oncology/Medical Imaging
- Clinical Relevance: Breast cancer metastasis detection

The task requires identifying whether the center 32x32px region of a patch 
contains at least one pixel of tumor tissue.
"""

from .grade import grade
from .prepare import prepare

__all__ = ["grade", "prepare"] 