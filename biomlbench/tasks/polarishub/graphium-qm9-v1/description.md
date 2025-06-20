# qm9-v1

## Background
QM9 is a well-known dataset in the field of 3D GNNs. It consists of 19 graph-level quantum properties associated to an energy-minimized 3D conformation of the molecules. It is considered a simple dataset since all the molecules have at most 9 heavy atoms. We chose QM9 in our ToyMix since it is very similar to the larger proposed quantum datasets, PCQM4M_multitask and PM6_83M, but with smaller molecules.

## Assay information
Computed geometric, energetic, electronic, and thermodynamic properties for 134k stable small organic molecules made up of CHONF. All properties were calculated at the B3LYP/6-31G(2df,p) level of quantum chemistry. For the predominant stoichiometry, C7H10O2, there are 6,095 constitutional isomers among the 134k molecules, with reported energies, enthalpies, and free energies of atomization at the more accurate G4MP2 level of theory.

## Benchmarking
**The goal** of this benchmark is to have the best predictive model for atomic coordinates and calculated properties.

## Data resource
Reference: [Quantum chemistry structures and properties of 134 kilo molecules](https://www.nature.com/articles/sdata201422)

---

**Source:** [Polaris Hub - graphium/qm9-v1](https://polarishub.io)  
**Task Type:** multi_task  
**Main Metric:** mean_squared_error

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `smiles`
- Target column: `{'u298_atom', 'h298_atom', 'cv', 'B', 'g298_atom', 'h298', 'alpha', 'gap', 'zpve', 'u298', 'mu', 'A', 'u0', 'u0_atom', 'g298', 'homo', 'lumo', 'r2', 'C'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `graphium/qm9-v1`.
Main metric: **mean_squared_error**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
