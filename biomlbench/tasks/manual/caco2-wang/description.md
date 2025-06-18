# Caco-2 Cell Permeability Prediction

## Background

The human colon epithelial cancer cell line, Caco-2, is used as an in vitro model to simulate the human intestinal tissue. The experimental result on the rate of drug passing through the Caco-2 cells can approximate the rate at which the drug permeates through the human intestinal tissue.

This is a critical property in drug discovery as it predicts how well a drug candidate will be absorbed in the gastrointestinal tract, which directly impacts bioavailability and therapeutic efficacy.

## Task Description

**Task Type**: Regression
**Input**: SMILES string representing molecular structure
**Output**: Caco-2 cell effective permeability (log Papp in cm/s)
**Metric**: Mean Absolute Error (lower is better)

Given a drug's SMILES string, predict its Caco-2 cell effective permeability. The permeability values are log-transformed apparent permeability coefficients (log Papp).

## Dataset Information

- **Training samples**: 728
- **Test samples**: 182
- **Feature**: SMILES (Simplified Molecular Input Line Entry System)
- **Target**: Caco-2 permeability (continuous values)

## Data Format

### Training Data (`train.csv`)
- `smiles`: SMILES string representation of the molecule
- `caco2_permeability`: Target permeability value (log Papp)

### Test Data (`test_features.csv`)
- `id`: Unique identifier for each test sample
- `smiles`: SMILES string representation of the molecule

### Submission Format (`submission.csv`)
- `id`: Must match the test sample IDs
- `caco2_permeability`: Your predicted permeability values

## Evaluation

Your predictions will be evaluated using Mean Absolute Error (MAE):

```
MAE = (1/n) * Σ|y_true - y_pred|
```

Lower scores are better. The current state-of-the-art approaches achieve MAE around 0.28-0.32.

## Baseline Approaches

Consider these molecular representation strategies:
1. **Molecular fingerprints**: ECFP, MACCS keys
2. **Molecular descriptors**: RDKit descriptors, Mordred
3. **Graph neural networks**: Operate directly on molecular graphs
4. **String-based models**: Transformer models on SMILES

## Clinical Relevance

Caco-2 permeability is crucial for:
- **Drug absorption prediction**: Estimates oral bioavailability
- **Lead optimization**: Guides medicinal chemistry decisions  
- **ADMET screening**: Early filtering of drug candidates
- **Regulatory submission**: Required data for drug approval

## References

Wang, NN et al. ADME Properties Evaluation in Drug Discovery: Prediction of Caco-2 Cell Permeability Using a Combination of NSGA-II and Boosting. Journal of Chemical Information and Modeling 2016, 56(4), 763-773.

Original benchmark: [Polaris Hub - tdcommons/caco2-wang](https://polarishub.io/benchmarks/tdcommons/caco2-wang) 