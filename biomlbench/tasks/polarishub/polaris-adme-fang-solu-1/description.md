# adme-fang-SOLU-1

![ADME](https://storage.googleapis.com/polaris-public/icons/icons8-whale-96-ADME.png) 

## Background

The goal of accessing ADME properties is to understand how a potential drug candidate interacts with the human body, including absorption, distribution, metabolism, and excretion. This knowledge is crucial for evaluating efficacy, safety, and clinical potential, guiding drug development for optimal therapeutic outcomes. Fang et al. 2023 has disclosed DMPK datasets collected over 20 months across six ADME in vitro endpoints, which are human and rat liver microsomal stability, MDR1-MDCK efflux ratio, solubility, and human and rat plasma protein binding. The dataset contains 885 to 3087 measures for the corresponding endpoints. The compounds show the chemical diversity across all ranges of the endpoints which are microsomal stability, plasma protein binding, permeability, and solubility.
## Benchmarking
**The goal** of this benchmark is to perform a single task, which is to the best predictive model for solubility. 


## Description of readout 
- **Readouts**: `LOG_SOLUBILITY`
- **Bioassay readout**: Solubility
- **Optimization objective**: Higher value
- **Number of data points**: train:  1578, test:  400


## Molecule data resource:
**Reference**: https://doi.org/10.1021/acs.jcim.3c00160

## Train/test split
In this benchmark set, the same train/test sets in the fang2023 paper were used for the 6 endpoints human and rat liver microsomal stability, MDR1-MDCK efflux ratio, solubility, and human and rat plasma protein binding, respectively. 
See more details at https://github.com/molecularinformatics/Computational-ADME/tree/main/MPNN.

**Distribution of the train/test in the chemical space**
![image](https://storage.googleapis.com/polaris-public/datasets/ADME/fang2023/figures/fang2023_ADME_public_v1_Sol_tsne_fang2023split.png)


## Related links
The full curation and creation process is documented [here](https://github.com/polaris-hub/polaris-recipes/blob/main/01_ADME).


---

**Source:** [Polaris Hub - polaris/adme-fang-solu-1](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** pearsonr

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `smiles`
- Target column: `{'LOG_SOLUBILITY'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `polaris/adme-fang-solu-1`.
Main metric: **pearsonr**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
