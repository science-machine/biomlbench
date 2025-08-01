"""
BioML-bench: A benchmark suite for evaluating machine learning agents on biomedical tasks.

This package provides a comprehensive framework for benchmarking AI agents on 
biomedical machine learning tasks including protein engineering, drug discovery, single cell omics,
medical imaging, and clinical biomarkers.
"""

__version__ = "0.1.0"
__author__ = "BioML-bench Team"

from biomlbench.registry import registry

__all__ = ["registry"] 