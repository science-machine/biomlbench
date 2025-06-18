# zinc12k-v1

## Background
ZINC is a free database of commercially-available compounds for virtual screening. ZINC contains over 230 million purchasable compounds in ready-to-dock, 3D formats. ZINC also contains over 750 million purchasable compounds that can be searched for analogs. ZINC12K contains a 12,000 sample subset of ZINC molecular graphs.

## Benchmarking
**The goal** of this benchmark is to have the best predictive model for synthetic accessibility (SA), logP and constrained solubility (Score).

## Description of readout:
- SA: Synthetic accessibility score.
- LogP: Log P, octanol-water partition coefficient.
- Score: constrained solubility which is the term logP − SA − cycle (octanol-water partition coefficients, logP, penalized by the synthetic accessibility score, SA, and number of long cycles, cycle).

The performance measure is the mean squared error (MSE) between the predicted and the ground truth value for each molecular graph.
  - Optimization objective: Lower value

## Data resource
References: 
- [ZINC: A Free Tool to Discover Chemistry for Biology](https://pubs.acs.org/doi/10.1021/ci3001277)
- [Benchmarking Graph Neural Networks](https://arxiv.org/pdf/2003.00982)

---

**Source:** [Polaris Hub - graphium/zinc12k-v1](https://polarishub.io)  
**Task Type:** multi_task  
**Main Metric:** mean_squared_error

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `smiles`
- Target column: `{'SA', 'logp', 'score'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `graphium/zinc12k-v1`.
Main metric: **mean_squared_error**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
