# adme-fang-r-1

![ADME](https://storage.googleapis.com/polaris-public/icons/icons8-whale-96-ADME.png) 

## Background

The goal of accessing ADME properties is to understand how a potential drug candidate interacts with the human body, including absorption, distribution, metabolism, and excretion. This knowledge is crucial for evaluating efficacy, safety, and clinical potential, guiding drug development for optimal therapeutic outcomes. [Fang et al. 2023](https://doi.org/10.1021/acs.jcim.3c00160) has disclosed DMPK datasets collected over 20 months across six ADME in vitro endpoints, which are human and rat liver microsomal stability, MDR1-MDCK efflux ratio, solubility, and human and rat plasma protein binding. The dataset contains 885 to 3087 measures for the corresponding endpoints. 

## Description of readout 
- **Readouts**: `LOG HLM_CLint (mL/min/kg)`, `LOG RLM_CLint (mL/min/kg)`, `LOG HPPB (mL/min/kg)`, `LOG RPPB (mL/min/kg)`, `LOG_MDR1-MDCK_ER`, `LOG_SOLUBILITY`
- **Bioassay readout**: Intrinsic clearance
- **Optimization objective**: Higher value
- **Number of data points**: train: 2812 test: 704

## Benchmarking
**The goal** of this benchmark is to perform a multitask learning, amd select the best models for predicting six adme endpoints altogether. 

## Molecule data resource:
**Reference**: https://doi.org/10.1021/acs.jcim.3c00160

## Train/test split
To discover more potential hits which are similar to the discovered hits, a random splitting was applied. 


**Distribution of the train/test in the chemical space**
![image](https://storage.googleapis.com/polaris-public/benchmarks/figures/fang2023_ADME_public_v1_tsne_random_split.png)

## Related links
The full curation and creation process is documented -> [notebook](https://github.com/polaris-hub/polaris-recipes/blob/main/02_MolProp/02_ADME_dataset.ipynb).

---

**Source:** [Polaris Hub - polaris/adme-fang-r-1](https://polarishub.io)  
**Task Type:** multi_task  
**Main Metric:** mean_squared_error

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `smiles`
- Target column: `{'LOG_RPPB', 'LOG_RLM_CLint', 'LOG_SOLUBILITY', 'LOG_HPPB', 'LOG_HLM_CLint', 'LOG_MDR1-MDCK_ER'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `polaris/adme-fang-r-1`.
Main metric: **mean_squared_error**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
