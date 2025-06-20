# bend-zeroshot-variant-effects-disease

# 🧬  [BEND](https://github.com/frederikkemarin/BEND) Zero-shot prediction of disease variants from ClinVar

Predicting variant effects is a binary problem, where single-bp mutations are classified as either having an effect or not. Each variant is a genomic position with a mutation $x ∈ {A, C, G, T}$ and a label $y ∈ {0, 1}$ indicating whether it is pathogenic or benign. The adjacent 512 bp serve as context. 

As this is a zero-shot task, we used the cosine distance in embedding space between a variant nucleotide and its reference nucleotide as the prediction score in BEND.


---

**Source:** [Polaris Hub - mlls/bend-zeroshot-variant-effects-disease](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** roc_auc

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `sequence_wt`
- Target column: `{'label'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `mlls/bend-zeroshot-variant-effects-disease`.
Main metric: **roc_auc**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
