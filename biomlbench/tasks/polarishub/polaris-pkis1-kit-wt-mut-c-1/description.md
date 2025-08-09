# pkis1-kit-wt-mut-c-1


### Background
**KIT** (Proto-oncogene c-KIT) receptor plays a crucial role in regulating cell growth, differentiation, and survival. It's particularly important in the development of blood cells, melanocytes (the cells that produce melanin, the pigment responsible for skin, hair, and eye color), and certain cells in the gut. Mutations in the KIT gene can lead to uncontrolled cell growth and contribute to the development of certain types of cancer, including gastrointestinal stromal tumors (GISTs) and some types of leukemia

### Benchmarking

- **KIT wild type**:  In certain cancers, KIT signaling can be activated by other receptors or mutations upstream in the signaling pathway. Targeting these upstream factors can indirectly impact KIT signaling and downstream effects. An example of this is seen in some cases of acute myeloid leukemia (AML) where KIT is expressed without mutations, but other upstream mutations can lead to aberrant KIT activation.
- **KIT selectivity**: This dataset includes KIT wild type and reported mutants `KIT T6701`, `KIT V560G`. D816V results in constitutive phosphorylation of Kit, activation of Stat5 signaling (PMID: 19865100, PMID: 18390729), induces mastocytosis and tumor formation in mice (PMID: 21148330) and confers resistance to Kit inhibitors (PMID: 22301675, PMID: 19164557). 

**The goal** of this benchmark is to perform a multitask, which is to the best predictive model for 
- Selectivity towards the mutants.
- Optimization of the bioactivity % inhibition.
- Discovery of potential hits in new chemical space.


## Description of readout:
- **Readouts**: `CLASS_KIT`, `CLASS_KIT_(T670I_mutant)`, `CLASS_KIT_(V560G_mutant)`
- **Bioassay readout**: percentage of inhnibition.
- **Number of data points**: train:  277 test:  87
- **Thresholds**: 
    - KIT_(T6701_mutant): > 75
    - KIT_(V560G_mutant): > 75
    - KIT: > 75
- **Optimization objective**: postive label (1)

### Data resource: 
- **Reference**: [PKIS1](https://pubmed.ncbi.nlm.nih.gov/26501955)

### Train/test split
Given the benchmarking goal, a scaffold-based splitting approach was applied to ensure training and test sets contain distinct chemical structures while maintaining the diversity of scaffolds.

**Distribution of the train/test in the chemical space**
![image](https://storage.googleapis.com/polaris-public/datasets/kinases/kit/figures/drewry_kit_wt_t670i_v560g_v1_tsne_scaffold_split.png)

**For more details of this benchmark** -> [notebook](https://github.com/polaris-hub/polaris-recipes/blob/main/03_Kinases/KIT)

## Related benchmarks
- polaris/drewry_kit_wt_t670i_v560g_clf_v1
> Note: It's recommanded to evaluate your methods agaisnt all the benchmarks related to this dataset. 


---

**Source:** [Polaris Hub - polaris/pkis1-kit-wt-mut-c-1](https://polarishub.io)  
**Task Type:** multi_task  
**Main Metric:** pr_auc

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `smiles`
- Target column: `{'CLASS_KIT_(T6701_mutant)', 'CLASS_KIT_(V560G_mutant)', 'CLASS_KIT'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `polaris/pkis1-kit-wt-mut-c-1`.
Main metric: **pr_auc**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
