# pkis2-egfr-wt-reg-v2

![molprop](https://storage.googleapis.com/polaris-public/icons/icons8-fox-60-kinases.png)


## Background
**EGFR (Epidermal Growth Factor Receptor) kinase** is a type of receptor tyrosine kinase that plays a significant role in cell growth, proliferation, and survival. Mutations or overexpression of EGFR have been associated with various diseases, particularly cancer.

## Benchmarking
 **EGFR Wild type**:  Targeting wild-type EGFR with small-molecule inhibitors, such as erlotinib, is an ongoing area of research in the treatment of glioblastoma. While early findings are promising, the complexity of glioblastoma biology presents challenges that require further investigation to improve treatment outcomes for patients.

 **EGFR L858R:** While EGFR TKIs initially demonstrate impressive responses in NSCLC patients with the L858R mutation, resistance to these drugs can develop over time. However, newer generations of EGFR TKIs, like osimertinib, have been developed to target these resistant mutations.

**The goal** of this benchmark is to select the best predictive model for 
- Optimization of the bioactivity % inhibition.
- Discovery of potential hits in new chemical space.

## Description of readout 
- **Readouts**: `EGFR`, `EGFR_L858R`
- **Bioassay readout**: Percentage of inhnibition.
- **Optimization objective**: Higher the more potent.
- **Number of data points**: train:  290 test:  74

## Data resource: 
- **Reference**: [PKIS1](https://pubmed.ncbi.nlm.nih.gov/26501955)


## Train/test split
Given the benchmarking goal, a scaffold-based splitting approach was applied to ensure training and test sets contain distinct chemical structures while maintaining the diversity of scaffolds.

**Distribution of the train/test in the chemical space**
![image](https://storage.googleapis.com/polaris-public/datasets/kinases/egfr/figures/egfr_wt_l858r_v1_tnse_scaffold_split.png)


## Related links
The full curation and creation process is documented -> [notebook](https://github.com/polaris-hub/polaris-recipes/blob/main/03_Kinases/EGFR)

## Related benchmarks
- polaris/drewry_egfr_wildtype_singletask_reg_v1
- polaris/drewry_egfr_wildtype_singletask_clf_v1
> Note: It's recommanded to evaluate your methods agaisnt all the benchmarks related to this dataset. 


---

**Source:** [Polaris Hub - polaris/pkis2-egfr-wt-reg-v2](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** mean_squared_error

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `MOL_smiles`
- Target column: `{'EGFR'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `polaris/pkis2-egfr-wt-reg-v2`.
Main metric: **mean_squared_error**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
