# pkis2-ret-wt-r-1

![molprop](https://storage.googleapis.com/polaris-public/icons/icons8-fox-60-kinases.png)

### Background
**RET (Rearranged during Transfection)** is a proto-oncogene that codes for a receptor tyrosine kinase. This means it produces a protein that plays a role in signaling pathways within cells, particularly related to cell growth and differentiation. When activated, RET helps regulate cell survival, proliferation, and differentiation. Mutations or alterations in the RET gene can lead to uncontrolled cell growth and potentially the development of cancer.

### Benchmarking
- **RET wild type**: In some cases, targeting both mutant and wild-type RET together can be more effective than targeting only one form as Combination Therapies. In certain cancer types, such as some subtypes of non-small cell lung cancer (NSCLC), the RET signaling pathway can interact with other oncogenic pathways, such as the EGFR (epidermal growth factor receptor) pathway. Targeting both pathways simultaneously might offer a synergistic effect and improve treatment outcomes.

**The goal** of this benchmark is to perform a single task, which is to the best predictive model for 
- Optimization of the bioactivity % inhibition for RET wile type.
- Discovery of potential hits in new chemical space.



### Description of readout 
- **Readouts**: `RET`
- **Bioassay readout**: percentage of inhnibition.
- **Optimization objective**: Higher value
- **Number of data points**: train:  534 test:  106


### Data resource: 
- **Reference**: [PKIS2](https://www.ncbi.nlm.nih.gov/pubmed/28767711)

### Train/test split
Given the benchmarking goal, a scaffold-based splitting approach was applied to ensure training and test sets contain distinct chemical structures while maintaining the diversity of scaffolds.

**Distribution of the train/test in the chemical space**
![image](https://storage.googleapis.com/polaris-public/datasets/kinases/ret/figures/drewry_ret_wildtype_v1_tnse_scaffold_split.png)

## Related links
The full curation and creation process is documented -> [notebook](https://github.com/polaris-hub/polaris-recipes/blob/main/03_Kinases/RET)

## Related benchmarks
- polaris/drewry_ret_wildtype_singletask_clf_v1
- polaris/drewry_ret_wt_v804l_y791f_multitask_reg_v1
- polaris/drewry_ret_wt_v804l_y791f_multitask_clf_v1
> Note: It's recommanded to evaluate your methods agaisnt all the benchmarks related to this dataset. 


---

**Source:** [Polaris Hub - polaris/pkis2-ret-wt-r-1](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** mean_squared_error

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `smiles`
- Target column: `{'RET'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `polaris/pkis2-ret-wt-r-1`.
Main metric: **mean_squared_error**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
