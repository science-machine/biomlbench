# EGFR_binders_binary_cls-v0

## Background

### Target details
Epidermal Growth Factor Receptor (EGFR) is a transmembrane protein that plays a critical role in cell growth, differentiation, and survival. It is frequently overexpressed or mutated in various cancers, including non-small cell lung cancer, colorectal cancer, and head and neck cancer. This makes EGFR a crucial target for cancer therapies such as Cetuximab, an antibody with more than 1B USD in annual revenue. 

- Target Protein: EGFR
- Organism: HUMAN
- Uniprot Accession ID: [P00533](https://www.uniprot.org/uniprotkb/P00533/entry)
- Protein sequence: LEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLENLQIIRGNMYYENSYALAVLSNYDANKTGLKELPMRNLQEILHGAVRFSNNPALCNVESIQWRDIVSSDFLSNMSMDFQNHLGSCQKCDPSCPNGSCWGAGEENCQKLTKIICAQQCSGRCRGKSPSDCCHNQCAAGCTGPRESDCLVCRKFRDEATCKDTCPPLMLYNPTTYQMDVNPEGKYSFGATCVKKCPRNYVVTDHGSCVRACGADSYEMEEDGVRKCKKCEGPCRKVCNGIGIGEFKDSLSINATNIKHFKNCTSISGDLHILPVAFRGDSFTHTPPLDPQELDILKTVKEITGFLLIQAWPENRTDLHAFENLEIIRGRTKQHGQFSLAVVSLNITSLGLRSLKEISDGDVIISGNKNLCYANTINWKKLFGTSGQKTKIISNRGENSCKATGQVCHALCSPEGCWGPEPRDCVSCRNVSRGRECVDKCNLLEGEPREFVENSECIQCHPECLPQAMNITCTGRGPDNCIQCAHYIDGPHCVKTCPAGVMGENNTLVWKYADAGHVCHLCHPNCTYGCTGPGLEGCPTNGPKIPS
- Structure PDB: [6ARU](https://www.rcsb.org/structure/6aru)


![64ru](https://cdn.rcsb.org/images/structures/6aru_assembly-1.jpeg)

### Binding protein designs
This dataset contains 202 designed EGFR-binding protein sequences, along with experimental binding affinity results tested by the AdaptyvBio team.

## Benchmark description

This retrospective benchmark evaluates protein design methods by challenging participants to design a binding protein for the extracellular domain of EGFR, a cancer-associated drug target. A set of 202 previously designed protein sequences, along with their experimental binding affinities (binary labels), is available for testing. 

`Balanced accuracy` is used to evaluate the performance of design methods in differentiating between binders and non-binders.

## Additional notes
- `Cetuximab_scFv` was used as the positive control of the binding assay.
- `ahmedsameh-Q3` and `ahmedsameh-yy2` were disqualified from the competition due to these similarity levels to the known EGFR-binder sequence. They are kpet in this benchmark solely for evaluation purposes.

## Reference: 
- https://design.adaptyvbio.com/
- https://foundry.adaptyvbio.com/egfr_design_competition
- https://github.com/adaptyvbio/egfr_competition_1

---

**Source:** [Polaris Hub - adaptyv-bio/egfr-binders-binary-cls-v0](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** balanced_accuracy

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `sequence`
- Target column: `{'binding_class'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `adaptyv-bio/egfr-binders-binary-cls-v0`.
Main metric: **balanced_accuracy**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
