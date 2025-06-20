# molprop-250k-reg-v2


![molprop](https://storage.googleapis.com/polaris-public/icons/icons8-bear-100-Molprop.png)

## Background

Molecular representations are crucial for understanding molecular structure, predicting properties, QSAR studies, toxicology and chemical modeling and other aspects in drug discovery tasks. Therefore, benchmarks for molecular representations are critical tools that drive progress in the field of computational chemistry and drug design. In recent years, many large models have been trained for learning molecular representation. The aim is to evaluate large pretrained models are capable of predicting various “easy-to-compute” molecular properties. 


## Benchmarking 

The objective is to comprehend the proficiency of a model in predicting these 'easy' properties, gauging its effectiveness. Ideally, any pre-trained models should, at the very least, demonstrate good performance in those tasks before applying them to the downstream tasks. 

## Description of readout
 The computed properties are molecular weight, fraction of sp3 carbon atoms (fsp3), number of rotatable bonds, topological polar surface area, computed logP, formal charge, number of charged atoms, refractivity and number of aromatic rings. These properties are widely used in molecule design and molecule prioritization.
 **Number of data points**: train:  202887, test:  46568

## Data resource
**Reference**: https://pubs.acs.org/doi/10.1021/acs.jcim.5b00559 

**Raw data**: https://raw.githubusercontent.com/aspuru-guzik-group/chemical_vae/master/models/zinc_properties/250k_randm_zinc_drugs_clean_3.csv 

## Train/test split
The objective is to comprehend the proficiency of a model in predicting these 'easy' properties. In order to select the predictive models which is able to generalize to new chemical space, a scaffold split is used to generate trian/test sets. 

**Distribution of the train/test in the chemical space**
![image](https://storage.googleapis.com/polaris-public/datasets/molprop/figures/molprop250k_v1_tsne_scaffold_split.png)


## Related links
The full curation and creation process is documented [here](https://github.com/polaris-hub/polaris-recipes/blob/main/02_MolProp).

## Related benchmarks
- [molprop250kleadlike_multitask_reg_v1](https://polarishub.io/benchmarks/polaris/molprop250kleadlike-multitask-reg-v1)

---

**Source:** [Polaris Hub - polaris/molprop-250k-reg-v2](https://polarishub.io)  
**Task Type:** multi_task  
**Main Metric:** mean_squared_error

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `MOL_smiles`
- Target column: `{'n_aromatic_rings', 'tpsa', 'n_rotatable_bonds', 'clogp', 'mw', 'n_charged_atoms', 'fsp3', 'formal_charge', 'refractivity'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `polaris/molprop-250k-reg-v2`.
Main metric: **mean_squared_error**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
