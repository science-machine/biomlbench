# tox21-v1

## Background
Tox21 is a well-known dataset for researchers in machine learning for drug discovery. The data set provided by the Tox21 Data Challenge included approximately 12 000 compounds. It consists of a multi-label classification task with 12 labels, with most labels missing and a strong imbalance towards the negative class. Each subchallenge required the prediction of a different type of toxicity. The sub-challenges were split between two panels: Seven of the twelve sub-challenges dealt with Nuclear Receptor (NR) signaling pathways, the remaining five with the Stress Response (SR) pathways. Nuclear receptors are important components in cell communication and control, and are involved in development,
metabolism and proliferation. Toxicity can also cause cellular stress which in term can
lead to apoptosis. Therefore the Tox21 data also includes five tasks on various stress response indicators. We chose Tox21 in our ToyMix since it is very similar to the larger proposed bioassay dataset, PCBA_1328_1564k both in terms of sparsity and imbalance and to the L1000 datasets in terms of imbalance.

## Assay information
In qHTS (quantitative high throughput screens), many thousands of compounds are screened in a single experiment across a broad concentration range in order to generate concentration–response curves. The method identifies compounds with a wide range of activities with a much lower false-positive or false-negative rate. The resulting concentration–response curves can be classified to rapidly identify actives and inactives with a variety of potencies and efficacies, producing rich data sets that can be mined for reliable biological activities.

## Benchmarking
**The goal** of this benchmark is to have the best predictive model for classifying molecules as having activity against toxicity-related proteins: nuclear receptors and stress response indicators.

## Data resource
Reference: [Improving the human hazard characterization of chemicals: a Tox21 update.](https://europepmc.org/article/MED/23603828)




---

**Source:** [Polaris Hub - graphium/tox21-v1](https://polarishub.io)  
**Task Type:** multi_task  
**Main Metric:** f1

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `smiles`
- Target column: `{'NR-ER-LBD', 'NR-AhR', 'NR-AR', 'SR-p53', 'SR-HSE', 'NR-PPAR-gamma', 'NR-ER', 'SR-ARE', 'NR-Aromatase', 'SR-MMP', 'SR-ATAD5', 'NR-AR-LBD'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `graphium/tox21-v1`.
Main metric: **f1**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
