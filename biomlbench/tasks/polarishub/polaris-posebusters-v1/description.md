# posebusters_v1

![logo](https://posebusters.readthedocs.io/en/latest/_static/logo_square.png)
## Background
The PoseBusters dataset set is a new set of carefully-selected publicly-available crystal complexes from the PDB. It is a diverse set of recent high-quality protein–ligand complexes which contain drug-like molecules. It only contains complexes released since 2021 and therefore does not contain any complexes present in the PDBbind General Set v2020 used to train many of the methods. 


[Buttenschoen et al.](https://pubs.rsc.org/en/content/articlehtml/2024/sc/d3sc04185a) outlines the steps used to select the 308 unique proteins and 308 unique ligands in the PoseBusters dataset. The complexes were downloaded from the PDB as MMTF files, and PyMOL was used to remove solvents and all occurrences of the ligand of interest. The proteins were saved with their cofactors in PDB format, while the ligands were saved in SDF format.

## Benchmark description

This is a zero-shot benchmark that contains only a test set of 308 proteisn and ligands for evaluation.

Posebusters offers a series of ligand checkers, known as 'Posebuster Checkers,' to filter out undesired docked ligand conformers. It is recommended to apply these filters before uploading results to the Polaris Hub.

Only the extracted ligand from the docking output should be uploaded for evaluation, ensuring that the target protein coordinates have been aligned with the original crystal structure.

This benchmark uses the metric RMSD coverage (≤2Å), which calculates the percentage of molecules with an RMSD of less than 2Å compared to the reference. For more details, see the documentation on [polaris](https://github.com/polaris-hub/polaris/blob/6e402f9d58d80d0ffe0b499cfef69a0c28c0427c/polaris/evaluate/metrics/docking_metrics.py#L39). .


## Data source
- Orignial: https://zenodo.org/records/8278563
- Polaris: gs://polaris-public/polaris-recipes/org-polaris/posebusters/posebusters_paper_data/posebusters_benchmark_set

## Other links
- Paper: https://pubs.rsc.org/en/content/articlelanding/2024/sc/d3sc04185a
- Github: https://github.com/maabuu/posebusters/tree/main

---

**Source:** [Polaris Hub - polaris/posebusters-v1](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** rmsd_coverage

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `protein`
- Target column: `{'ligand'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `polaris/posebusters-v1`.
Main metric: **rmsd_coverage**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
