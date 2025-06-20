# adme-fang-reg-v1

![ADME](https://storage.googleapis.com/polaris-public/icons/icons8-whale-96-ADME.png) 

## Background

The goal of assessing ADME properties is to understand how a potential drug candidate interacts with the human body, including absorption, distribution, metabolism, and excretion. This understanding is essential for evaluating the drug's efficacy, safety, and clinical potential, guiding drug development towards optimal therapeutic outcomes. [Fang et al. (2023)](https://doi.org/10.1021/acs.jcim.3c00160) disclosed DMPK datasets collected over 20 months, covering six ADME in vitro endpoints: human and rat liver microsomal stability, MDR1-MDCK efflux ratio, solubility, and human and rat plasma protein binding. The dataset includes between 885 and 3,087 measurements for each corresponding endpoint.


## Description of readout 
- **Readouts**: `LOG HLM_CLint (mL/min/kg)`, `LOG RLM_CLint (mL/min/kg)`, `LOG HPPB (mL/min/kg)`, `LOG RPPB (mL/min/kg)`, `LOG_MDR1-MDCK_ER`, `LOG_SOLUBILITY`
- **Bioassay readout**: Intrinsic clearance
- **Optimization objective**: Higher value

## Benchmarking
**The goal** of this benchmark is to perform multitask learning and select the best models for predicting six ADME endpoints altogether.

## Train/test split
To discover more potential hits that are similar to the discovered hits, a random splitting was applied.

**Distribution of the train/test in the chemical space**
![image](https://storage.googleapis.com/polaris-public/biogen/fang2023_ADME/figures/multitask_random_chemspace.png)


---

**Source:** [Polaris Hub - biogen/adme-fang-reg-v1](https://polarishub.io)  
**Task Type:** multi_task  
**Main Metric:** mean_squared_error

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `MOL_smiles`
- Target column: `{'LOG_RPPB', 'LOG_RLM_CLint', 'LOG_SOLUBILITY', 'LOG_HPPB', 'LOG_HLM_CLint', 'LOG_MDR1-MDCK_ER'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `biogen/adme-fang-reg-v1`.
Main metric: **mean_squared_error**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
