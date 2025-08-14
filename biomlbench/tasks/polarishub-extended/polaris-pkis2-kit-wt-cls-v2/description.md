# pkis2-kit-wt-cls-v2


### Background
**KIT** (Proto-oncogene c-KIT) receptor plays a crucial role in regulating cell growth, differentiation, and survival. It's particularly important in the development of blood cells, melanocytes (the cells that produce melanin, the pigment responsible for skin, hair, and eye color), and certain cells in the gut. Mutations in the KIT gene can lead to uncontrolled cell growth and contribute to the development of certain types of cancer, including gastrointestinal stromal tumors (GISTs) and some types of leukemia

### Benchmarking

- **KIT wild type**:  In certain cancers, KIT signaling can be activated by other receptors or mutations upstream in the signaling pathway. Targeting these upstream factors can indirectly impact KIT signaling and downstream effects. An example of this is seen in some cases of acute myeloid leukemia (AML) where KIT is expressed without mutations, but other upstream mutations can lead to aberrant KIT activation.

**The goal** of this benchmark is to perform a multitask, which is to the best predictive model for 
- Selectivity towards the mutants.
- Optimization of the bioactivity % inhibition.
- Discovery of potential hits in new chemical space.


## Description of readout:
- **Readouts**: `CLS_KIT`
- **Bioassay readout**: percentage of inhnibition.
- **Thresholds**: 
    - KIT: > 80
- **Optimization objective**: postive label (1)

### Data resource: 
- **Reference**: [PKIS1](https://pubmed.ncbi.nlm.nih.gov/26501955)

### Train/test split
Given the benchmarking goal, a scaffold-based splitting approach was applied to ensure training and test sets contain distinct chemical structures while maintaining the diversity of scaffolds.

**Distribution of the train/test in the chemical space**
![image](https://storage.googleapis.com/polaris-public/polaris-recipes/org-polaris/drewry2017_pkis2_subset/figures/scaffold_split_chemspace.png)

**For more details of this benchmark** -> [notebook](https://github.com/polaris-hub/polaris-recipes/blob/main/org-Polaris/drewry2017_pkis2_subset/benchmarks/02_pkis2-kit_wt_benchmark.ipynb)



---

**Source:** [Polaris Hub - polaris/pkis2-kit-wt-cls-v2](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** pr_auc

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `MOL_smiles`
- Target column: `{'CLS_KIT'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `polaris/pkis2-kit-wt-cls-v2`.
Main metric: **pr_auc**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
