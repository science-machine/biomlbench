# adme-fang-RPPB-1

**Benchmarking goal:** 

Single tasks for the six endpoints: As the original paper, author established regression tasks for each ADME endpoints with predefined train-test set for the model training. In this benchmark set, the same train/test sets in the fang2023 paper were used for the 6 endpoints human and rat liver microsomal stability, MDR1-MDCK efflux ratio, solubility, and human and rat plasma protein binding, respectively. 

## Benchmarking
*The goal** of this benchmark is to perform a single task, which is to the best predictive model for human liver microsomal stability. 


## Description of readout 
- **Readouts**: `LOG PLASMA PROTEIN BINDING (RAT) (% unbound)`
- **Bioassay readout**: Rat plasma protein binding 
- **Optimization objective**: Lower value


## Molecule data resource:
**Reference**: https://doi.org/10.1021/acs.jcim.3c00160

## Train/test split
In this benchmark set, the same train/test sets in the fang2023 paper were used for the 6 endpoints human and rat liver microsomal stability, MDR1-MDCK efflux ratio, solubility, and human and rat plasma protein binding, respectively. 
See more details at https://github.com/molecularinformatics/Computational-ADME/tree/main/MPNN.

**Distribution of the train/test in the chemical space**
![image](https://storage.googleapis.com/polaris-public/datasets/ADME/fang2023/figures/fang2023_ADME_public_v1_hPPB_tsne_fang2023split.png)

## Related links
The full curation and creation process is documented [here](https://github.com/polaris-hub/polaris-recipes/blob/main/01_ADME).


---

**Source:** [Polaris Hub - polaris/adme-fang-rppb-1](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** pearsonr

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `smiles`
- Target column: `{'LOG_RPPB'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `polaris/adme-fang-rppb-1`.
Main metric: **pearsonr**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
