# CLI Reference

BioML-bench's command-line interface for managing biomedical benchmark tasks and running agent evaluations.

::: biomlbench.cli
    options:
      show_source: false
      heading_level: 2

## Main Help

```
usage: cli.py [-h] {prepare,grade,grade-sample,dev,run-baseline,run-agent} ...

Runs agents on biomedical ML tasks.

positional arguments:
  {prepare,grade,grade-sample,dev,run-baseline,run-agent}
                        Sub-command to run.
    prepare             Download and prepare tasks for the BioML-bench
                        dataset.
    grade               Grade a submission to the eval, comprising of several
                        task submissions
    grade-sample        Grade a single sample (dataset) in the eval
    dev                 Developer tools for extending BioML-bench.
    run-baseline        Run a baseline agent on a biomedical task
    run-agent           Run an AI agent on biomedical tasks

options:
  -h, --help            show this help message and exit

```

## Commands

### `biomlbench prepare`

```
usage: cli.py prepare [-h] [-t TASK_ID] [-d DATASET_ID] [-a] [--lite]
                      [-l LIST] [--domain DOMAIN] [--task-type TASK_TYPE]
                      [--keep-raw] [--data-dir DATA_DIR]
                      [--overwrite-checksums] [--overwrite-leaderboard]
                      [--skip-verification]

options:
  -h, --help            show this help message and exit
  -t TASK_ID, --task-id TASK_ID
                        ID of the task to prepare in 'folder/task' format.
                        Examples: manual/caco2-wang, polarishub/tdcommons-
                        admet
  -d DATASET_ID, --dataset-id DATASET_ID
                        ID of the dataset to prepare in 'folder/dataset'
                        format. Examples: proteingym-
                        dms/A0A247D711_LISMN_Stadelmann_2021
  -a, --all             Prepare all tasks.
  --lite                Prepare all the low difficulty tasks (BioML-bench
                        Lite).
  -l LIST, --list LIST  Prepare a list of tasks specified line by line in a
                        text file.
  --domain DOMAIN       Prepare all tasks for a specific biomedical domain
                        (e.g., oncology, drug_discovery).
  --task-type TASK_TYPE
                        Prepare all tasks of a specific type (e.g.,
                        medical_imaging, protein_engineering).
  --keep-raw            Keep the raw task files after the task has been
                        prepared.
  --data-dir DATA_DIR   Path to the directory where the data will be stored.
  --overwrite-checksums
                        [For Developers] Overwrite the checksums file for the
                        task.
  --overwrite-leaderboard
                        [For Developers] Overwrite the leaderboard file for
                        the task.
  --skip-verification   [For Developers] Skip the verification of the
                        checksums.

```

### `biomlbench run-agent`

```
usage: cli.py run-agent [-h] --agent AGENT [--task-id TASK_ID]
                        [--task-list TASK_LIST] [--n-workers N_WORKERS]
                        [--n-seeds N_SEEDS]
                        [--container-config CONTAINER_CONFIG]
                        [--retain-container] [--data-dir DATA_DIR]
                        [--output-dir OUTPUT_DIR]

options:
  -h, --help            show this help message and exit
  --agent AGENT         Agent ID to run (e.g., dummy, aide, aide/dev)
  --task-id TASK_ID     Single task ID to run. Valid options:
                        ['manual/caco2-wang',
                        'openproblems/batch_integration', 'polarishub/adaptyv-
                        bio-egfr-binders-binary-cls-v0', 'polarishub/biogen-
                        adme-fang-hclint-reg-v1', 'polarishub/biogen-adme-
                        fang-hppb-reg-v1', 'polarishub/biogen-adme-fang-perm-
                        reg-v1', 'polarishub/biogen-adme-fang-rclint-reg-v1',
                        'polarishub/biogen-adme-fang-reg-v1',
                        'polarishub/biogen-adme-fang-rppb-reg-v1',
                        'polarishub/biogen-adme-fang-solu-reg-v1',
                        'polarishub/graphium-l1000-mcf7-v1',
                        'polarishub/graphium-l1000-vcap-v1',
                        'polarishub/graphium-pcba-1328-1564k-v1',
                        'polarishub/graphium-qm9-v1', 'polarishub/graphium-
                        tox21-v1', 'polarishub/graphium-zinc12k-v1',
                        'polarishub/mlls-bend-zeroshot-variant-effects-
                        disease', 'polarishub/mlls-bend-zeroshot-variant-
                        effects-expression', 'polarishub/molecularml-
                        moleculeace-chembl1862-ki', 'polarishub/molecularml-
                        moleculeace-chembl1871-ki', 'polarishub/molecularml-
                        moleculeace-chembl2034-ki', 'polarishub/molecularml-
                        moleculeace-chembl204-ki', 'polarishub/molecularml-
                        moleculeace-chembl2047-ec50', 'polarishub/molecularml-
                        moleculeace-chembl214-ki', 'polarishub/molecularml-
                        moleculeace-chembl2147-ki', 'polarishub/molecularml-
                        moleculeace-chembl218-ec50', 'polarishub/molecularml-
                        moleculeace-chembl219-ki', 'polarishub/molecularml-
                        moleculeace-chembl228-ki', 'polarishub/molecularml-
                        moleculeace-chembl231-ki', 'polarishub/molecularml-
                        moleculeace-chembl233-ki', 'polarishub/molecularml-
                        moleculeace-chembl234-ki', 'polarishub/molecularml-
                        moleculeace-chembl235-ec50', 'polarishub/molecularml-
                        moleculeace-chembl236-ki', 'polarishub/molecularml-
                        moleculeace-chembl237-ec50', 'polarishub/molecularml-
                        moleculeace-chembl237-ki', 'polarishub/molecularml-
                        moleculeace-chembl238-ki', 'polarishub/molecularml-
                        moleculeace-chembl239-ec50', 'polarishub/molecularml-
                        moleculeace-chembl244-ki', 'polarishub/molecularml-
                        moleculeace-chembl262-ki', 'polarishub/molecularml-
                        moleculeace-chembl264-ki', 'polarishub/molecularml-
                        moleculeace-chembl2835-ki', 'polarishub/molecularml-
                        moleculeace-chembl287-ki', 'polarishub/molecularml-
                        moleculeace-chembl2971-ki', 'polarishub/molecularml-
                        moleculeace-chembl3979-ec50', 'polarishub/molecularml-
                        moleculeace-chembl4005-ki', 'polarishub/molecularml-
                        moleculeace-chembl4203-ki', 'polarishub/molecularml-
                        moleculeace-chembl4616-ec50', 'polarishub/molecularml-
                        moleculeace-chembl4792-ki', 'polarishub/novartis-adme-
                        novartis-cyp3a4-cls', 'polarishub/novartis-adme-
                        novartis-cyp3a4-reg', 'polarishub/polaris-adme-fang-
                        hclint-1', 'polarishub/polaris-adme-fang-hppb-1',
                        'polarishub/polaris-adme-fang-perm-1',
                        'polarishub/polaris-adme-fang-r-1',
                        'polarishub/polaris-adme-fang-rclint-1',
                        'polarishub/polaris-adme-fang-rppb-1',
                        'polarishub/polaris-adme-fang-solu-1',
                        'polarishub/polaris-hello-world-benchmark',
                        'polarishub/polaris-molprop-250k-leadlike-r-1',
                        'polarishub/polaris-molprop-250k-r-1',
                        'polarishub/polaris-molprop-250k-reg-v2',
                        'polarishub/polaris-molprop-leadlike-250k-reg-v2',
                        'polarishub/polaris-pkis1-egfr-wt-mut-c-1',
                        'polarishub/polaris-pkis1-egfr-wt-mut-r-1',
                        'polarishub/polaris-pkis1-kit-wt-mut-c-1',
                        'polarishub/polaris-pkis1-kit-wt-mut-r-1',
                        'polarishub/polaris-pkis1-ret-wt-mut-c-1',
                        'polarishub/polaris-pkis1-ret-wt-mut-r-1',
                        'polarishub/polaris-pkis2-egfr-wt-c-1',
                        'polarishub/polaris-pkis2-egfr-wt-r-1',
                        'polarishub/polaris-pkis2-egfr-wt-reg-v2',
                        'polarishub/polaris-pkis2-kit-wt-c-1',
                        'polarishub/polaris-pkis2-kit-wt-cls-v2',
                        'polarishub/polaris-pkis2-kit-wt-r-1',
                        'polarishub/polaris-pkis2-kit-wt-reg-v2',
                        'polarishub/polaris-pkis2-lok-slk-c-1',
                        'polarishub/polaris-pkis2-lok-slk-cls-v2',
                        'polarishub/polaris-pkis2-lok-slk-r-1',
                        'polarishub/polaris-pkis2-ret-wt-c-1',
                        'polarishub/polaris-pkis2-ret-wt-cls-v2',
                        'polarishub/polaris-pkis2-ret-wt-r-1',
                        'polarishub/polaris-pkis2-ret-wt-reg-v2',
                        'polarishub/polaris-posebusters-v1',
                        'polarishub/recursion-rxrx-compound-gene-activity-
                        benchmark', 'polarishub/tdcommons-ames',
                        'polarishub/tdcommons-bbb-martins',
                        'polarishub/tdcommons-bioavailability-ma',
                        'polarishub/tdcommons-caco2-wang',
                        'polarishub/tdcommons-clearance-hepatocyte-az',
                        'polarishub/tdcommons-clearance-microsome-az',
                        'polarishub/tdcommons-cyp2c9-substrate-carbonmangels',
                        'polarishub/tdcommons-cyp2c9-veith',
                        'polarishub/tdcommons-cyp2d6-substrate-carbonmangels',
                        'polarishub/tdcommons-cyp2d6-veith',
                        'polarishub/tdcommons-cyp3a4-substrate-carbonmangels',
                        'polarishub/tdcommons-cyp3a4-veith',
                        'polarishub/tdcommons-dili', 'polarishub/tdcommons-
                        half-life-obach', 'polarishub/tdcommons-herg',
                        'polarishub/tdcommons-hia-hou', 'polarishub/tdcommons-
                        ld50-zhu', 'polarishub/tdcommons-lipophilicity-
                        astrazeneca', 'polarishub/tdcommons-pgp-broccatelli',
                        'polarishub/tdcommons-ppbr-az', 'polarishub/tdcommons-
                        solubility-aqsoldb', 'polarishub/tdcommons-vdss-
                        lombardo', 'polarishub/vishrut64-cho-dna-expression-
                        prediction-dataset-task', 'polarishub/vishrut64-rna-
                        expression-prediction-dataset-task', 'proteingym-
                        dms/A0A140D2T1_ZIKV_Sourisseau_2019', 'proteingym-
                        dms/A0A192B1T2_9HIV1_Haddox_2018', 'proteingym-
                        dms/A0A1I9GEU1_NEIME_Kennouche_2019', 'proteingym-
                        dms/A0A247D711_LISMN_Stadelmann_2021', 'proteingym-
                        dms/A0A2Z5U3Z0_9INFA_Doud_2016', 'proteingym-
                        dms/A0A2Z5U3Z0_9INFA_Wu_2014', 'proteingym-
                        dms/A4D664_9INFA_Soh_2019', 'proteingym-
                        dms/A4GRB6_PSEAI_Chen_2020', 'proteingym-
                        dms/AACC1_PSEAI_Dandage_2018', 'proteingym-
                        dms/ACE2_HUMAN_Chan_2020', 'proteingym-
                        dms/ADRB2_HUMAN_Jones_2020', 'proteingym-
                        dms/AICDA_HUMAN_Gajula_2014_3cycles', 'proteingym-
                        dms/AMIE_PSEAE_Wrenbeck_2017', 'proteingym-
                        dms/ANCSZ_Hobbs_2022', 'proteingym-
                        dms/ARGR_ECOLI_Tsuboyama_2023_1AOY', 'proteingym-
                        dms/B2L11_HUMAN_Dutta_2010_binding-Mcl-1',
                        'proteingym-dms/BLAT_ECOLX_Deng_2012', 'proteingym-
                        dms/BLAT_ECOLX_Firnberg_2014', 'proteingym-
                        dms/BLAT_ECOLX_Jacquier_2013', 'proteingym-
                        dms/BLAT_ECOLX_Stiffler_2015', 'proteingym-
                        dms/BRCA1_HUMAN_Findlay_2018', 'proteingym-
                        dms/BRCA2_HUMAN_Erwood_2022_HEK293T', 'proteingym-
                        dms/C6KNH7_9INFA_Lee_2018', 'proteingym-
                        dms/CALM1_HUMAN_Weile_2017', 'proteingym-
                        dms/CAR11_HUMAN_Meitlis_2020_gof', 'proteingym-
                        dms/CAR11_HUMAN_Meitlis_2020_lof', 'proteingym-
                        dms/CAS9_STRP1_Spencer_2017_positive', 'proteingym-
                        dms/CASP3_HUMAN_Roychowdhury_2020', 'proteingym-
                        dms/CASP7_HUMAN_Roychowdhury_2020', 'proteingym-
                        dms/CBS_HUMAN_Sun_2020', 'proteingym-
                        dms/CCDB_ECOLI_Adkar_2012', 'proteingym-
                        dms/CCDB_ECOLI_Tripathi_2016', 'proteingym-
                        dms/CCR5_HUMAN_Gill_2023', 'proteingym-
                        dms/CD19_HUMAN_Klesmith_2019_FMC_singles',
                        'proteingym-dms/CP2C9_HUMAN_Amorosi_2021_abundance',
                        'proteingym-dms/CP2C9_HUMAN_Amorosi_2021_activity',
                        'proteingym-dms/DLG4_RAT_McLaughlin_2012',
                        'proteingym-dms/DN7A_SACS2_Tsuboyama_2023_1JIC',
                        'proteingym-dms/DYR_ECOLI_Nguyen_2023', 'proteingym-
                        dms/DYR_ECOLI_Thompson_2019', 'proteingym-
                        dms/ENVZ_ECOLI_Ghose_2023', 'proteingym-
                        dms/ENV_HV1B9_DuenasDecamp_2016', 'proteingym-
                        dms/ENV_HV1BR_Haddox_2016', 'proteingym-
                        dms/ERBB2_HUMAN_Elazar_2016', 'proteingym-
                        dms/ESTA_BACSU_Nutschel_2020', 'proteingym-
                        dms/FKBP3_HUMAN_Tsuboyama_2023_2KFV', 'proteingym-
                        dms/GAL4_YEAST_Kitzman_2015', 'proteingym-
                        dms/GDIA_HUMAN_Silverstein_2021', 'proteingym-
                        dms/GLPA_HUMAN_Elazar_2016', 'proteingym-
                        dms/HCP_LAMBD_Tsuboyama_2023_2L6Q', 'proteingym-
                        dms/HEM3_HUMAN_Loggerenberg_2023', 'proteingym-
                        dms/HMDH_HUMAN_Jiang_2019', 'proteingym-
                        dms/HSP82_YEAST_Cote-Hammarlof_2020_growth-H2O2',
                        'proteingym-dms/HSP82_YEAST_Flynn_2019', 'proteingym-
                        dms/HSP82_YEAST_Mishra_2016', 'proteingym-
                        dms/HXK4_HUMAN_Gersing_2022_activity', 'proteingym-
                        dms/HXK4_HUMAN_Gersing_2023_abundance', 'proteingym-
                        dms/I6TAH8_I68A0_Doud_2015', 'proteingym-
                        dms/IF1_ECOLI_Kelsic_2016', 'proteingym-
                        dms/ILF3_HUMAN_Tsuboyama_2023_2L33', 'proteingym-
                        dms/KCNE1_HUMAN_Muhammad_2023_expression',
                        'proteingym-dms/KCNE1_HUMAN_Muhammad_2023_function',
                        'proteingym-dms/KCNH2_HUMAN_Kozek_2020', 'proteingym-
                        dms/KCNJ2_MOUSE_Coyote-Maestas_2022_function',
                        'proteingym-dms/KCNJ2_MOUSE_Coyote-
                        Maestas_2022_surface', 'proteingym-
                        dms/KKA2_KLEPN_Melnikov_2014', 'proteingym-
                        dms/LGK_LIPST_Klesmith_2015', 'proteingym-
                        dms/LYAM1_HUMAN_Elazar_2016', 'proteingym-
                        dms/MET_HUMAN_Estevam_2023', 'proteingym-
                        dms/MK01_HUMAN_Brenan_2016', 'proteingym-
                        dms/MLAC_ECOLI_MacRae_2023', 'proteingym-
                        dms/MSH2_HUMAN_Jia_2020', 'proteingym-
                        dms/MTH3_HAEAE_RockahShmuel_2015', 'proteingym-
                        dms/MTHR_HUMAN_Weile_2021', 'proteingym-
                        dms/NCAP_I34A1_Doud_2015', 'proteingym-
                        dms/NPC1_HUMAN_Erwood_2022_HEK293T', 'proteingym-
                        dms/NPC1_HUMAN_Erwood_2022_RPE1', 'proteingym-
                        dms/NRAM_I33A0_Jiang_2016', 'proteingym-
                        dms/NUD15_HUMAN_Suiter_2020', 'proteingym-
                        dms/OPSD_HUMAN_Wan_2019', 'proteingym-
                        dms/OTC_HUMAN_Lo_2023', 'proteingym-
                        dms/OTU7A_HUMAN_Tsuboyama_2023_2L2D', 'proteingym-
                        dms/OXDA_RHOTO_Vanella_2023_activity', 'proteingym-
                        dms/OXDA_RHOTO_Vanella_2023_expression', 'proteingym-
                        dms/P53_HUMAN_Giacomelli_2018_Null_Etoposide',
                        'proteingym-
                        dms/P53_HUMAN_Giacomelli_2018_Null_Nutlin',
                        'proteingym-dms/P53_HUMAN_Giacomelli_2018_WT_Nutlin',
                        'proteingym-dms/P53_HUMAN_Kotler_2018', 'proteingym-
                        dms/P84126_THETH_Chan_2017', 'proteingym-
                        dms/PAI1_HUMAN_Huttinger_2021', 'proteingym-
                        dms/PA_I34A1_Wu_2015', 'proteingym-
                        dms/PKN1_HUMAN_Tsuboyama_2023_1URF', 'proteingym-
                        dms/POLG_CXB3N_Mattenberger_2021', 'proteingym-
                        dms/POLG_DEN26_Suphatrakul_2023', 'proteingym-
                        dms/POLG_HCVJF_Qi_2014', 'proteingym-
                        dms/PPARG_HUMAN_Majithia_2016', 'proteingym-
                        dms/PPM1D_HUMAN_Miller_2022', 'proteingym-
                        dms/PRKN_HUMAN_Clausen_2023', 'proteingym-
                        dms/PTEN_HUMAN_Matreyek_2021', 'proteingym-
                        dms/PTEN_HUMAN_Mighell_2018', 'proteingym-
                        dms/Q2N0S5_9HIV1_Haddox_2018', 'proteingym-
                        dms/Q53Z42_HUMAN_McShan_2019_binding-TAPBPR',
                        'proteingym-dms/Q53Z42_HUMAN_McShan_2019_expression',
                        'proteingym-dms/Q59976_STRSQ_Romero_2015',
                        'proteingym-dms/Q837P4_ENTFA_Meier_2023', 'proteingym-
                        dms/Q837P5_ENTFA_Meier_2023', 'proteingym-
                        dms/R1AB_SARS2_Flynn_2022', 'proteingym-
                        dms/RAF1_HUMAN_Zinkus-Boltz_2019', 'proteingym-
                        dms/RASH_HUMAN_Bandaru_2017', 'proteingym-
                        dms/RDRP_I33A0_Li_2023', 'proteingym-
                        dms/REV_HV1H2_Fernandes_2016', 'proteingym-
                        dms/RL40A_YEAST_Mavor_2016', 'proteingym-
                        dms/RL40A_YEAST_Roscoe_2013', 'proteingym-
                        dms/RL40A_YEAST_Roscoe_2014', 'proteingym-
                        dms/RNC_ECOLI_Weeks_2023', 'proteingym-
                        dms/RPC1_LAMBD_Li_2019_high-expression', 'proteingym-
                        dms/RPC1_LAMBD_Li_2019_low-expression', 'proteingym-
                        dms/RS15_GEOSE_Tsuboyama_2023_1A32', 'proteingym-
                        dms/S22A1_HUMAN_Yee_2023_abundance', 'proteingym-
                        dms/S22A1_HUMAN_Yee_2023_activity', 'proteingym-
                        dms/SBI_STAAM_Tsuboyama_2023_2JVG', 'proteingym-
                        dms/SC6A4_HUMAN_Young_2021', 'proteingym-
                        dms/SCIN_STAAR_Tsuboyama_2023_2QFF', 'proteingym-
                        dms/SCN5A_HUMAN_Glazer_2019', 'proteingym-
                        dms/SERC_HUMAN_Xie_2023', 'proteingym-
                        dms/SHOC2_HUMAN_Kwon_2022', 'proteingym-
                        dms/SOX30_HUMAN_Tsuboyama_2023_7JJK', 'proteingym-
                        dms/SPIKE_SARS2_Starr_2020_binding', 'proteingym-
                        dms/SPIKE_SARS2_Starr_2020_expression', 'proteingym-
                        dms/SQSTM_MOUSE_Tsuboyama_2023_2RRU', 'proteingym-
                        dms/SRC_HUMAN_Ahler_2019', 'proteingym-
                        dms/SRC_HUMAN_Chakraborty_2023_binding-DAS_25uM',
                        'proteingym-dms/SRC_HUMAN_Nguyen_2022', 'proteingym-
                        dms/SUMO1_HUMAN_Weile_2017', 'proteingym-
                        dms/SYUA_HUMAN_Newberry_2020', 'proteingym-
                        dms/TADBP_HUMAN_Bolognesi_2019', 'proteingym-
                        dms/TAT_HV1BR_Fernandes_2016', 'proteingym-
                        dms/TPK1_HUMAN_Weile_2017', 'proteingym-
                        dms/TPMT_HUMAN_Matreyek_2018', 'proteingym-
                        dms/TPOR_HUMAN_Bridgford_2020', 'proteingym-
                        dms/TRPC_SACS2_Chan_2017', 'proteingym-
                        dms/TRPC_THEMA_Chan_2017', 'proteingym-
                        dms/UBC9_HUMAN_Weile_2017', 'proteingym-
                        dms/UBE4B_MOUSE_Starita_2013', 'proteingym-
                        dms/VG08_BPP22_Tsuboyama_2023_2GP8', 'proteingym-
                        dms/VKOR1_HUMAN_Chiasson_2020_abundance', 'proteingym-
                        dms/VKOR1_HUMAN_Chiasson_2020_activity', 'proteingym-
                        dms/VRPI_BPT7_Tsuboyama_2023_2WNM']
  --task-list TASK_LIST
                        Path to text file with task IDs (one per line) for
                        multi-task runs
  --n-workers N_WORKERS
                        Number of parallel workers for multi-task runs
  --n-seeds N_SEEDS     Number of random seeds to run per task
  --container-config CONTAINER_CONFIG
                        Path to JSON file with Docker container configuration
  --retain-container    Keep container after run for debugging
  --data-dir DATA_DIR   Path to the directory where task data is stored
  --output-dir OUTPUT_DIR
                        Directory to save agent run outputs

```

### `biomlbench grade`

```
usage: cli.py grade [-h] --submission SUBMISSION --output-dir OUTPUT_DIR
                    [--data-dir DATA_DIR]

options:
  -h, --help            show this help message and exit
  --submission SUBMISSION
                        Path to the JSONL file of submissions. Refer to
                        README.md#submission-format for the required format.
  --output-dir OUTPUT_DIR
                        Path to the directory where the evaluation metrics
                        will be saved.
  --data-dir DATA_DIR   Path to the directory where the data used for grading
                        is stored.

```

### `biomlbench grade-sample`

```
usage: cli.py grade-sample [-h] [--data-dir DATA_DIR]
                           submission task_id dataset_id

positional arguments:
  submission           Path to the submission file.
  task_id              ID of the task to grade in 'folder/task' format.
                       Examples: manual/caco2-wang
  dataset_id           ID of the dataset to grade. Should be a folder in the
                       task's directory.

options:
  -h, --help           show this help message and exit
  --data-dir DATA_DIR  Path to the directory where the data will be stored.

```

### `biomlbench run-baseline`

```
usage: cli.py run-baseline [-h] [--baseline BASELINE]
                           [--output-dir OUTPUT_DIR] [--data-dir DATA_DIR]
                           [--seed SEED]
                           task_id

positional arguments:
  task_id               ID of the task to run baseline on. Valid options:
                        ['manual/caco2-wang',
                        'openproblems/batch_integration', 'polarishub/adaptyv-
                        bio-egfr-binders-binary-cls-v0', 'polarishub/biogen-
                        adme-fang-hclint-reg-v1', 'polarishub/biogen-adme-
                        fang-hppb-reg-v1', 'polarishub/biogen-adme-fang-perm-
                        reg-v1', 'polarishub/biogen-adme-fang-rclint-reg-v1',
                        'polarishub/biogen-adme-fang-reg-v1',
                        'polarishub/biogen-adme-fang-rppb-reg-v1',
                        'polarishub/biogen-adme-fang-solu-reg-v1',
                        'polarishub/graphium-l1000-mcf7-v1',
                        'polarishub/graphium-l1000-vcap-v1',
                        'polarishub/graphium-pcba-1328-1564k-v1',
                        'polarishub/graphium-qm9-v1', 'polarishub/graphium-
                        tox21-v1', 'polarishub/graphium-zinc12k-v1',
                        'polarishub/mlls-bend-zeroshot-variant-effects-
                        disease', 'polarishub/mlls-bend-zeroshot-variant-
                        effects-expression', 'polarishub/molecularml-
                        moleculeace-chembl1862-ki', 'polarishub/molecularml-
                        moleculeace-chembl1871-ki', 'polarishub/molecularml-
                        moleculeace-chembl2034-ki', 'polarishub/molecularml-
                        moleculeace-chembl204-ki', 'polarishub/molecularml-
                        moleculeace-chembl2047-ec50', 'polarishub/molecularml-
                        moleculeace-chembl214-ki', 'polarishub/molecularml-
                        moleculeace-chembl2147-ki', 'polarishub/molecularml-
                        moleculeace-chembl218-ec50', 'polarishub/molecularml-
                        moleculeace-chembl219-ki', 'polarishub/molecularml-
                        moleculeace-chembl228-ki', 'polarishub/molecularml-
                        moleculeace-chembl231-ki', 'polarishub/molecularml-
                        moleculeace-chembl233-ki', 'polarishub/molecularml-
                        moleculeace-chembl234-ki', 'polarishub/molecularml-
                        moleculeace-chembl235-ec50', 'polarishub/molecularml-
                        moleculeace-chembl236-ki', 'polarishub/molecularml-
                        moleculeace-chembl237-ec50', 'polarishub/molecularml-
                        moleculeace-chembl237-ki', 'polarishub/molecularml-
                        moleculeace-chembl238-ki', 'polarishub/molecularml-
                        moleculeace-chembl239-ec50', 'polarishub/molecularml-
                        moleculeace-chembl244-ki', 'polarishub/molecularml-
                        moleculeace-chembl262-ki', 'polarishub/molecularml-
                        moleculeace-chembl264-ki', 'polarishub/molecularml-
                        moleculeace-chembl2835-ki', 'polarishub/molecularml-
                        moleculeace-chembl287-ki', 'polarishub/molecularml-
                        moleculeace-chembl2971-ki', 'polarishub/molecularml-
                        moleculeace-chembl3979-ec50', 'polarishub/molecularml-
                        moleculeace-chembl4005-ki', 'polarishub/molecularml-
                        moleculeace-chembl4203-ki', 'polarishub/molecularml-
                        moleculeace-chembl4616-ec50', 'polarishub/molecularml-
                        moleculeace-chembl4792-ki', 'polarishub/novartis-adme-
                        novartis-cyp3a4-cls', 'polarishub/novartis-adme-
                        novartis-cyp3a4-reg', 'polarishub/polaris-adme-fang-
                        hclint-1', 'polarishub/polaris-adme-fang-hppb-1',
                        'polarishub/polaris-adme-fang-perm-1',
                        'polarishub/polaris-adme-fang-r-1',
                        'polarishub/polaris-adme-fang-rclint-1',
                        'polarishub/polaris-adme-fang-rppb-1',
                        'polarishub/polaris-adme-fang-solu-1',
                        'polarishub/polaris-hello-world-benchmark',
                        'polarishub/polaris-molprop-250k-leadlike-r-1',
                        'polarishub/polaris-molprop-250k-r-1',
                        'polarishub/polaris-molprop-250k-reg-v2',
                        'polarishub/polaris-molprop-leadlike-250k-reg-v2',
                        'polarishub/polaris-pkis1-egfr-wt-mut-c-1',
                        'polarishub/polaris-pkis1-egfr-wt-mut-r-1',
                        'polarishub/polaris-pkis1-kit-wt-mut-c-1',
                        'polarishub/polaris-pkis1-kit-wt-mut-r-1',
                        'polarishub/polaris-pkis1-ret-wt-mut-c-1',
                        'polarishub/polaris-pkis1-ret-wt-mut-r-1',
                        'polarishub/polaris-pkis2-egfr-wt-c-1',
                        'polarishub/polaris-pkis2-egfr-wt-r-1',
                        'polarishub/polaris-pkis2-egfr-wt-reg-v2',
                        'polarishub/polaris-pkis2-kit-wt-c-1',
                        'polarishub/polaris-pkis2-kit-wt-cls-v2',
                        'polarishub/polaris-pkis2-kit-wt-r-1',
                        'polarishub/polaris-pkis2-kit-wt-reg-v2',
                        'polarishub/polaris-pkis2-lok-slk-c-1',
                        'polarishub/polaris-pkis2-lok-slk-cls-v2',
                        'polarishub/polaris-pkis2-lok-slk-r-1',
                        'polarishub/polaris-pkis2-ret-wt-c-1',
                        'polarishub/polaris-pkis2-ret-wt-cls-v2',
                        'polarishub/polaris-pkis2-ret-wt-r-1',
                        'polarishub/polaris-pkis2-ret-wt-reg-v2',
                        'polarishub/polaris-posebusters-v1',
                        'polarishub/recursion-rxrx-compound-gene-activity-
                        benchmark', 'polarishub/tdcommons-ames',
                        'polarishub/tdcommons-bbb-martins',
                        'polarishub/tdcommons-bioavailability-ma',
                        'polarishub/tdcommons-caco2-wang',
                        'polarishub/tdcommons-clearance-hepatocyte-az',
                        'polarishub/tdcommons-clearance-microsome-az',
                        'polarishub/tdcommons-cyp2c9-substrate-carbonmangels',
                        'polarishub/tdcommons-cyp2c9-veith',
                        'polarishub/tdcommons-cyp2d6-substrate-carbonmangels',
                        'polarishub/tdcommons-cyp2d6-veith',
                        'polarishub/tdcommons-cyp3a4-substrate-carbonmangels',
                        'polarishub/tdcommons-cyp3a4-veith',
                        'polarishub/tdcommons-dili', 'polarishub/tdcommons-
                        half-life-obach', 'polarishub/tdcommons-herg',
                        'polarishub/tdcommons-hia-hou', 'polarishub/tdcommons-
                        ld50-zhu', 'polarishub/tdcommons-lipophilicity-
                        astrazeneca', 'polarishub/tdcommons-pgp-broccatelli',
                        'polarishub/tdcommons-ppbr-az', 'polarishub/tdcommons-
                        solubility-aqsoldb', 'polarishub/tdcommons-vdss-
                        lombardo', 'polarishub/vishrut64-cho-dna-expression-
                        prediction-dataset-task', 'polarishub/vishrut64-rna-
                        expression-prediction-dataset-task', 'proteingym-
                        dms/A0A140D2T1_ZIKV_Sourisseau_2019', 'proteingym-
                        dms/A0A192B1T2_9HIV1_Haddox_2018', 'proteingym-
                        dms/A0A1I9GEU1_NEIME_Kennouche_2019', 'proteingym-
                        dms/A0A247D711_LISMN_Stadelmann_2021', 'proteingym-
                        dms/A0A2Z5U3Z0_9INFA_Doud_2016', 'proteingym-
                        dms/A0A2Z5U3Z0_9INFA_Wu_2014', 'proteingym-
                        dms/A4D664_9INFA_Soh_2019', 'proteingym-
                        dms/A4GRB6_PSEAI_Chen_2020', 'proteingym-
                        dms/AACC1_PSEAI_Dandage_2018', 'proteingym-
                        dms/ACE2_HUMAN_Chan_2020', 'proteingym-
                        dms/ADRB2_HUMAN_Jones_2020', 'proteingym-
                        dms/AICDA_HUMAN_Gajula_2014_3cycles', 'proteingym-
                        dms/AMIE_PSEAE_Wrenbeck_2017', 'proteingym-
                        dms/ANCSZ_Hobbs_2022', 'proteingym-
                        dms/ARGR_ECOLI_Tsuboyama_2023_1AOY', 'proteingym-
                        dms/B2L11_HUMAN_Dutta_2010_binding-Mcl-1',
                        'proteingym-dms/BLAT_ECOLX_Deng_2012', 'proteingym-
                        dms/BLAT_ECOLX_Firnberg_2014', 'proteingym-
                        dms/BLAT_ECOLX_Jacquier_2013', 'proteingym-
                        dms/BLAT_ECOLX_Stiffler_2015', 'proteingym-
                        dms/BRCA1_HUMAN_Findlay_2018', 'proteingym-
                        dms/BRCA2_HUMAN_Erwood_2022_HEK293T', 'proteingym-
                        dms/C6KNH7_9INFA_Lee_2018', 'proteingym-
                        dms/CALM1_HUMAN_Weile_2017', 'proteingym-
                        dms/CAR11_HUMAN_Meitlis_2020_gof', 'proteingym-
                        dms/CAR11_HUMAN_Meitlis_2020_lof', 'proteingym-
                        dms/CAS9_STRP1_Spencer_2017_positive', 'proteingym-
                        dms/CASP3_HUMAN_Roychowdhury_2020', 'proteingym-
                        dms/CASP7_HUMAN_Roychowdhury_2020', 'proteingym-
                        dms/CBS_HUMAN_Sun_2020', 'proteingym-
                        dms/CCDB_ECOLI_Adkar_2012', 'proteingym-
                        dms/CCDB_ECOLI_Tripathi_2016', 'proteingym-
                        dms/CCR5_HUMAN_Gill_2023', 'proteingym-
                        dms/CD19_HUMAN_Klesmith_2019_FMC_singles',
                        'proteingym-dms/CP2C9_HUMAN_Amorosi_2021_abundance',
                        'proteingym-dms/CP2C9_HUMAN_Amorosi_2021_activity',
                        'proteingym-dms/DLG4_RAT_McLaughlin_2012',
                        'proteingym-dms/DN7A_SACS2_Tsuboyama_2023_1JIC',
                        'proteingym-dms/DYR_ECOLI_Nguyen_2023', 'proteingym-
                        dms/DYR_ECOLI_Thompson_2019', 'proteingym-
                        dms/ENVZ_ECOLI_Ghose_2023', 'proteingym-
                        dms/ENV_HV1B9_DuenasDecamp_2016', 'proteingym-
                        dms/ENV_HV1BR_Haddox_2016', 'proteingym-
                        dms/ERBB2_HUMAN_Elazar_2016', 'proteingym-
                        dms/ESTA_BACSU_Nutschel_2020', 'proteingym-
                        dms/FKBP3_HUMAN_Tsuboyama_2023_2KFV', 'proteingym-
                        dms/GAL4_YEAST_Kitzman_2015', 'proteingym-
                        dms/GDIA_HUMAN_Silverstein_2021', 'proteingym-
                        dms/GLPA_HUMAN_Elazar_2016', 'proteingym-
                        dms/HCP_LAMBD_Tsuboyama_2023_2L6Q', 'proteingym-
                        dms/HEM3_HUMAN_Loggerenberg_2023', 'proteingym-
                        dms/HMDH_HUMAN_Jiang_2019', 'proteingym-
                        dms/HSP82_YEAST_Cote-Hammarlof_2020_growth-H2O2',
                        'proteingym-dms/HSP82_YEAST_Flynn_2019', 'proteingym-
                        dms/HSP82_YEAST_Mishra_2016', 'proteingym-
                        dms/HXK4_HUMAN_Gersing_2022_activity', 'proteingym-
                        dms/HXK4_HUMAN_Gersing_2023_abundance', 'proteingym-
                        dms/I6TAH8_I68A0_Doud_2015', 'proteingym-
                        dms/IF1_ECOLI_Kelsic_2016', 'proteingym-
                        dms/ILF3_HUMAN_Tsuboyama_2023_2L33', 'proteingym-
                        dms/KCNE1_HUMAN_Muhammad_2023_expression',
                        'proteingym-dms/KCNE1_HUMAN_Muhammad_2023_function',
                        'proteingym-dms/KCNH2_HUMAN_Kozek_2020', 'proteingym-
                        dms/KCNJ2_MOUSE_Coyote-Maestas_2022_function',
                        'proteingym-dms/KCNJ2_MOUSE_Coyote-
                        Maestas_2022_surface', 'proteingym-
                        dms/KKA2_KLEPN_Melnikov_2014', 'proteingym-
                        dms/LGK_LIPST_Klesmith_2015', 'proteingym-
                        dms/LYAM1_HUMAN_Elazar_2016', 'proteingym-
                        dms/MET_HUMAN_Estevam_2023', 'proteingym-
                        dms/MK01_HUMAN_Brenan_2016', 'proteingym-
                        dms/MLAC_ECOLI_MacRae_2023', 'proteingym-
                        dms/MSH2_HUMAN_Jia_2020', 'proteingym-
                        dms/MTH3_HAEAE_RockahShmuel_2015', 'proteingym-
                        dms/MTHR_HUMAN_Weile_2021', 'proteingym-
                        dms/NCAP_I34A1_Doud_2015', 'proteingym-
                        dms/NPC1_HUMAN_Erwood_2022_HEK293T', 'proteingym-
                        dms/NPC1_HUMAN_Erwood_2022_RPE1', 'proteingym-
                        dms/NRAM_I33A0_Jiang_2016', 'proteingym-
                        dms/NUD15_HUMAN_Suiter_2020', 'proteingym-
                        dms/OPSD_HUMAN_Wan_2019', 'proteingym-
                        dms/OTC_HUMAN_Lo_2023', 'proteingym-
                        dms/OTU7A_HUMAN_Tsuboyama_2023_2L2D', 'proteingym-
                        dms/OXDA_RHOTO_Vanella_2023_activity', 'proteingym-
                        dms/OXDA_RHOTO_Vanella_2023_expression', 'proteingym-
                        dms/P53_HUMAN_Giacomelli_2018_Null_Etoposide',
                        'proteingym-
                        dms/P53_HUMAN_Giacomelli_2018_Null_Nutlin',
                        'proteingym-dms/P53_HUMAN_Giacomelli_2018_WT_Nutlin',
                        'proteingym-dms/P53_HUMAN_Kotler_2018', 'proteingym-
                        dms/P84126_THETH_Chan_2017', 'proteingym-
                        dms/PAI1_HUMAN_Huttinger_2021', 'proteingym-
                        dms/PA_I34A1_Wu_2015', 'proteingym-
                        dms/PKN1_HUMAN_Tsuboyama_2023_1URF', 'proteingym-
                        dms/POLG_CXB3N_Mattenberger_2021', 'proteingym-
                        dms/POLG_DEN26_Suphatrakul_2023', 'proteingym-
                        dms/POLG_HCVJF_Qi_2014', 'proteingym-
                        dms/PPARG_HUMAN_Majithia_2016', 'proteingym-
                        dms/PPM1D_HUMAN_Miller_2022', 'proteingym-
                        dms/PRKN_HUMAN_Clausen_2023', 'proteingym-
                        dms/PTEN_HUMAN_Matreyek_2021', 'proteingym-
                        dms/PTEN_HUMAN_Mighell_2018', 'proteingym-
                        dms/Q2N0S5_9HIV1_Haddox_2018', 'proteingym-
                        dms/Q53Z42_HUMAN_McShan_2019_binding-TAPBPR',
                        'proteingym-dms/Q53Z42_HUMAN_McShan_2019_expression',
                        'proteingym-dms/Q59976_STRSQ_Romero_2015',
                        'proteingym-dms/Q837P4_ENTFA_Meier_2023', 'proteingym-
                        dms/Q837P5_ENTFA_Meier_2023', 'proteingym-
                        dms/R1AB_SARS2_Flynn_2022', 'proteingym-
                        dms/RAF1_HUMAN_Zinkus-Boltz_2019', 'proteingym-
                        dms/RASH_HUMAN_Bandaru_2017', 'proteingym-
                        dms/RDRP_I33A0_Li_2023', 'proteingym-
                        dms/REV_HV1H2_Fernandes_2016', 'proteingym-
                        dms/RL40A_YEAST_Mavor_2016', 'proteingym-
                        dms/RL40A_YEAST_Roscoe_2013', 'proteingym-
                        dms/RL40A_YEAST_Roscoe_2014', 'proteingym-
                        dms/RNC_ECOLI_Weeks_2023', 'proteingym-
                        dms/RPC1_LAMBD_Li_2019_high-expression', 'proteingym-
                        dms/RPC1_LAMBD_Li_2019_low-expression', 'proteingym-
                        dms/RS15_GEOSE_Tsuboyama_2023_1A32', 'proteingym-
                        dms/S22A1_HUMAN_Yee_2023_abundance', 'proteingym-
                        dms/S22A1_HUMAN_Yee_2023_activity', 'proteingym-
                        dms/SBI_STAAM_Tsuboyama_2023_2JVG', 'proteingym-
                        dms/SC6A4_HUMAN_Young_2021', 'proteingym-
                        dms/SCIN_STAAR_Tsuboyama_2023_2QFF', 'proteingym-
                        dms/SCN5A_HUMAN_Glazer_2019', 'proteingym-
                        dms/SERC_HUMAN_Xie_2023', 'proteingym-
                        dms/SHOC2_HUMAN_Kwon_2022', 'proteingym-
                        dms/SOX30_HUMAN_Tsuboyama_2023_7JJK', 'proteingym-
                        dms/SPIKE_SARS2_Starr_2020_binding', 'proteingym-
                        dms/SPIKE_SARS2_Starr_2020_expression', 'proteingym-
                        dms/SQSTM_MOUSE_Tsuboyama_2023_2RRU', 'proteingym-
                        dms/SRC_HUMAN_Ahler_2019', 'proteingym-
                        dms/SRC_HUMAN_Chakraborty_2023_binding-DAS_25uM',
                        'proteingym-dms/SRC_HUMAN_Nguyen_2022', 'proteingym-
                        dms/SUMO1_HUMAN_Weile_2017', 'proteingym-
                        dms/SYUA_HUMAN_Newberry_2020', 'proteingym-
                        dms/TADBP_HUMAN_Bolognesi_2019', 'proteingym-
                        dms/TAT_HV1BR_Fernandes_2016', 'proteingym-
                        dms/TPK1_HUMAN_Weile_2017', 'proteingym-
                        dms/TPMT_HUMAN_Matreyek_2018', 'proteingym-
                        dms/TPOR_HUMAN_Bridgford_2020', 'proteingym-
                        dms/TRPC_SACS2_Chan_2017', 'proteingym-
                        dms/TRPC_THEMA_Chan_2017', 'proteingym-
                        dms/UBC9_HUMAN_Weile_2017', 'proteingym-
                        dms/UBE4B_MOUSE_Starita_2013', 'proteingym-
                        dms/VG08_BPP22_Tsuboyama_2023_2GP8', 'proteingym-
                        dms/VKOR1_HUMAN_Chiasson_2020_abundance', 'proteingym-
                        dms/VKOR1_HUMAN_Chiasson_2020_activity', 'proteingym-
                        dms/VRPI_BPT7_Tsuboyama_2023_2WNM']

options:
  -h, --help            show this help message and exit
  --baseline BASELINE   Baseline to run (simple, random, or task-specific
                        types like linear, rf, fingerprint). Use 'all' to run
                        all available baselines for the task.
  --output-dir OUTPUT_DIR
                        Directory to save baseline submissions
  --data-dir DATA_DIR   Path to the directory where the data is stored.
  --seed SEED           Random seed for reproducible baselines

```

### `biomlbench dev`

```
usage: cli.py dev [-h] {download-leaderboard,prepare-human-baselines} ...

positional arguments:
  {download-leaderboard,prepare-human-baselines}
                        Developer command to run.
    download-leaderboard
                        Download the leaderboard for a task.
    prepare-human-baselines
                        Prepare human baseline data for tasks.

options:
  -h, --help            show this help message and exit

```

#### `biomlbench dev download-leaderboard`

```
usage: cli.py dev download-leaderboard [-h] [-t TASK_ID] [--all] [--force]

options:
  -h, --help            show this help message and exit
  -t TASK_ID, --task-id TASK_ID
                        Name of the task to download the leaderboard for.
                        Valid options: ['manual/caco2-wang',
                        'openproblems/batch_integration', 'polarishub/adaptyv-
                        bio-egfr-binders-binary-cls-v0', 'polarishub/biogen-
                        adme-fang-hclint-reg-v1', 'polarishub/biogen-adme-
                        fang-hppb-reg-v1', 'polarishub/biogen-adme-fang-perm-
                        reg-v1', 'polarishub/biogen-adme-fang-rclint-reg-v1',
                        'polarishub/biogen-adme-fang-reg-v1',
                        'polarishub/biogen-adme-fang-rppb-reg-v1',
                        'polarishub/biogen-adme-fang-solu-reg-v1',
                        'polarishub/graphium-l1000-mcf7-v1',
                        'polarishub/graphium-l1000-vcap-v1',
                        'polarishub/graphium-pcba-1328-1564k-v1',
                        'polarishub/graphium-qm9-v1', 'polarishub/graphium-
                        tox21-v1', 'polarishub/graphium-zinc12k-v1',
                        'polarishub/mlls-bend-zeroshot-variant-effects-
                        disease', 'polarishub/mlls-bend-zeroshot-variant-
                        effects-expression', 'polarishub/molecularml-
                        moleculeace-chembl1862-ki', 'polarishub/molecularml-
                        moleculeace-chembl1871-ki', 'polarishub/molecularml-
                        moleculeace-chembl2034-ki', 'polarishub/molecularml-
                        moleculeace-chembl204-ki', 'polarishub/molecularml-
                        moleculeace-chembl2047-ec50', 'polarishub/molecularml-
                        moleculeace-chembl214-ki', 'polarishub/molecularml-
                        moleculeace-chembl2147-ki', 'polarishub/molecularml-
                        moleculeace-chembl218-ec50', 'polarishub/molecularml-
                        moleculeace-chembl219-ki', 'polarishub/molecularml-
                        moleculeace-chembl228-ki', 'polarishub/molecularml-
                        moleculeace-chembl231-ki', 'polarishub/molecularml-
                        moleculeace-chembl233-ki', 'polarishub/molecularml-
                        moleculeace-chembl234-ki', 'polarishub/molecularml-
                        moleculeace-chembl235-ec50', 'polarishub/molecularml-
                        moleculeace-chembl236-ki', 'polarishub/molecularml-
                        moleculeace-chembl237-ec50', 'polarishub/molecularml-
                        moleculeace-chembl237-ki', 'polarishub/molecularml-
                        moleculeace-chembl238-ki', 'polarishub/molecularml-
                        moleculeace-chembl239-ec50', 'polarishub/molecularml-
                        moleculeace-chembl244-ki', 'polarishub/molecularml-
                        moleculeace-chembl262-ki', 'polarishub/molecularml-
                        moleculeace-chembl264-ki', 'polarishub/molecularml-
                        moleculeace-chembl2835-ki', 'polarishub/molecularml-
                        moleculeace-chembl287-ki', 'polarishub/molecularml-
                        moleculeace-chembl2971-ki', 'polarishub/molecularml-
                        moleculeace-chembl3979-ec50', 'polarishub/molecularml-
                        moleculeace-chembl4005-ki', 'polarishub/molecularml-
                        moleculeace-chembl4203-ki', 'polarishub/molecularml-
                        moleculeace-chembl4616-ec50', 'polarishub/molecularml-
                        moleculeace-chembl4792-ki', 'polarishub/novartis-adme-
                        novartis-cyp3a4-cls', 'polarishub/novartis-adme-
                        novartis-cyp3a4-reg', 'polarishub/polaris-adme-fang-
                        hclint-1', 'polarishub/polaris-adme-fang-hppb-1',
                        'polarishub/polaris-adme-fang-perm-1',
                        'polarishub/polaris-adme-fang-r-1',
                        'polarishub/polaris-adme-fang-rclint-1',
                        'polarishub/polaris-adme-fang-rppb-1',
                        'polarishub/polaris-adme-fang-solu-1',
                        'polarishub/polaris-hello-world-benchmark',
                        'polarishub/polaris-molprop-250k-leadlike-r-1',
                        'polarishub/polaris-molprop-250k-r-1',
                        'polarishub/polaris-molprop-250k-reg-v2',
                        'polarishub/polaris-molprop-leadlike-250k-reg-v2',
                        'polarishub/polaris-pkis1-egfr-wt-mut-c-1',
                        'polarishub/polaris-pkis1-egfr-wt-mut-r-1',
                        'polarishub/polaris-pkis1-kit-wt-mut-c-1',
                        'polarishub/polaris-pkis1-kit-wt-mut-r-1',
                        'polarishub/polaris-pkis1-ret-wt-mut-c-1',
                        'polarishub/polaris-pkis1-ret-wt-mut-r-1',
                        'polarishub/polaris-pkis2-egfr-wt-c-1',
                        'polarishub/polaris-pkis2-egfr-wt-r-1',
                        'polarishub/polaris-pkis2-egfr-wt-reg-v2',
                        'polarishub/polaris-pkis2-kit-wt-c-1',
                        'polarishub/polaris-pkis2-kit-wt-cls-v2',
                        'polarishub/polaris-pkis2-kit-wt-r-1',
                        'polarishub/polaris-pkis2-kit-wt-reg-v2',
                        'polarishub/polaris-pkis2-lok-slk-c-1',
                        'polarishub/polaris-pkis2-lok-slk-cls-v2',
                        'polarishub/polaris-pkis2-lok-slk-r-1',
                        'polarishub/polaris-pkis2-ret-wt-c-1',
                        'polarishub/polaris-pkis2-ret-wt-cls-v2',
                        'polarishub/polaris-pkis2-ret-wt-r-1',
                        'polarishub/polaris-pkis2-ret-wt-reg-v2',
                        'polarishub/polaris-posebusters-v1',
                        'polarishub/recursion-rxrx-compound-gene-activity-
                        benchmark', 'polarishub/tdcommons-ames',
                        'polarishub/tdcommons-bbb-martins',
                        'polarishub/tdcommons-bioavailability-ma',
                        'polarishub/tdcommons-caco2-wang',
                        'polarishub/tdcommons-clearance-hepatocyte-az',
                        'polarishub/tdcommons-clearance-microsome-az',
                        'polarishub/tdcommons-cyp2c9-substrate-carbonmangels',
                        'polarishub/tdcommons-cyp2c9-veith',
                        'polarishub/tdcommons-cyp2d6-substrate-carbonmangels',
                        'polarishub/tdcommons-cyp2d6-veith',
                        'polarishub/tdcommons-cyp3a4-substrate-carbonmangels',
                        'polarishub/tdcommons-cyp3a4-veith',
                        'polarishub/tdcommons-dili', 'polarishub/tdcommons-
                        half-life-obach', 'polarishub/tdcommons-herg',
                        'polarishub/tdcommons-hia-hou', 'polarishub/tdcommons-
                        ld50-zhu', 'polarishub/tdcommons-lipophilicity-
                        astrazeneca', 'polarishub/tdcommons-pgp-broccatelli',
                        'polarishub/tdcommons-ppbr-az', 'polarishub/tdcommons-
                        solubility-aqsoldb', 'polarishub/tdcommons-vdss-
                        lombardo', 'polarishub/vishrut64-cho-dna-expression-
                        prediction-dataset-task', 'polarishub/vishrut64-rna-
                        expression-prediction-dataset-task', 'proteingym-
                        dms/A0A140D2T1_ZIKV_Sourisseau_2019', 'proteingym-
                        dms/A0A192B1T2_9HIV1_Haddox_2018', 'proteingym-
                        dms/A0A1I9GEU1_NEIME_Kennouche_2019', 'proteingym-
                        dms/A0A247D711_LISMN_Stadelmann_2021', 'proteingym-
                        dms/A0A2Z5U3Z0_9INFA_Doud_2016', 'proteingym-
                        dms/A0A2Z5U3Z0_9INFA_Wu_2014', 'proteingym-
                        dms/A4D664_9INFA_Soh_2019', 'proteingym-
                        dms/A4GRB6_PSEAI_Chen_2020', 'proteingym-
                        dms/AACC1_PSEAI_Dandage_2018', 'proteingym-
                        dms/ACE2_HUMAN_Chan_2020', 'proteingym-
                        dms/ADRB2_HUMAN_Jones_2020', 'proteingym-
                        dms/AICDA_HUMAN_Gajula_2014_3cycles', 'proteingym-
                        dms/AMIE_PSEAE_Wrenbeck_2017', 'proteingym-
                        dms/ANCSZ_Hobbs_2022', 'proteingym-
                        dms/ARGR_ECOLI_Tsuboyama_2023_1AOY', 'proteingym-
                        dms/B2L11_HUMAN_Dutta_2010_binding-Mcl-1',
                        'proteingym-dms/BLAT_ECOLX_Deng_2012', 'proteingym-
                        dms/BLAT_ECOLX_Firnberg_2014', 'proteingym-
                        dms/BLAT_ECOLX_Jacquier_2013', 'proteingym-
                        dms/BLAT_ECOLX_Stiffler_2015', 'proteingym-
                        dms/BRCA1_HUMAN_Findlay_2018', 'proteingym-
                        dms/BRCA2_HUMAN_Erwood_2022_HEK293T', 'proteingym-
                        dms/C6KNH7_9INFA_Lee_2018', 'proteingym-
                        dms/CALM1_HUMAN_Weile_2017', 'proteingym-
                        dms/CAR11_HUMAN_Meitlis_2020_gof', 'proteingym-
                        dms/CAR11_HUMAN_Meitlis_2020_lof', 'proteingym-
                        dms/CAS9_STRP1_Spencer_2017_positive', 'proteingym-
                        dms/CASP3_HUMAN_Roychowdhury_2020', 'proteingym-
                        dms/CASP7_HUMAN_Roychowdhury_2020', 'proteingym-
                        dms/CBS_HUMAN_Sun_2020', 'proteingym-
                        dms/CCDB_ECOLI_Adkar_2012', 'proteingym-
                        dms/CCDB_ECOLI_Tripathi_2016', 'proteingym-
                        dms/CCR5_HUMAN_Gill_2023', 'proteingym-
                        dms/CD19_HUMAN_Klesmith_2019_FMC_singles',
                        'proteingym-dms/CP2C9_HUMAN_Amorosi_2021_abundance',
                        'proteingym-dms/CP2C9_HUMAN_Amorosi_2021_activity',
                        'proteingym-dms/DLG4_RAT_McLaughlin_2012',
                        'proteingym-dms/DN7A_SACS2_Tsuboyama_2023_1JIC',
                        'proteingym-dms/DYR_ECOLI_Nguyen_2023', 'proteingym-
                        dms/DYR_ECOLI_Thompson_2019', 'proteingym-
                        dms/ENVZ_ECOLI_Ghose_2023', 'proteingym-
                        dms/ENV_HV1B9_DuenasDecamp_2016', 'proteingym-
                        dms/ENV_HV1BR_Haddox_2016', 'proteingym-
                        dms/ERBB2_HUMAN_Elazar_2016', 'proteingym-
                        dms/ESTA_BACSU_Nutschel_2020', 'proteingym-
                        dms/FKBP3_HUMAN_Tsuboyama_2023_2KFV', 'proteingym-
                        dms/GAL4_YEAST_Kitzman_2015', 'proteingym-
                        dms/GDIA_HUMAN_Silverstein_2021', 'proteingym-
                        dms/GLPA_HUMAN_Elazar_2016', 'proteingym-
                        dms/HCP_LAMBD_Tsuboyama_2023_2L6Q', 'proteingym-
                        dms/HEM3_HUMAN_Loggerenberg_2023', 'proteingym-
                        dms/HMDH_HUMAN_Jiang_2019', 'proteingym-
                        dms/HSP82_YEAST_Cote-Hammarlof_2020_growth-H2O2',
                        'proteingym-dms/HSP82_YEAST_Flynn_2019', 'proteingym-
                        dms/HSP82_YEAST_Mishra_2016', 'proteingym-
                        dms/HXK4_HUMAN_Gersing_2022_activity', 'proteingym-
                        dms/HXK4_HUMAN_Gersing_2023_abundance', 'proteingym-
                        dms/I6TAH8_I68A0_Doud_2015', 'proteingym-
                        dms/IF1_ECOLI_Kelsic_2016', 'proteingym-
                        dms/ILF3_HUMAN_Tsuboyama_2023_2L33', 'proteingym-
                        dms/KCNE1_HUMAN_Muhammad_2023_expression',
                        'proteingym-dms/KCNE1_HUMAN_Muhammad_2023_function',
                        'proteingym-dms/KCNH2_HUMAN_Kozek_2020', 'proteingym-
                        dms/KCNJ2_MOUSE_Coyote-Maestas_2022_function',
                        'proteingym-dms/KCNJ2_MOUSE_Coyote-
                        Maestas_2022_surface', 'proteingym-
                        dms/KKA2_KLEPN_Melnikov_2014', 'proteingym-
                        dms/LGK_LIPST_Klesmith_2015', 'proteingym-
                        dms/LYAM1_HUMAN_Elazar_2016', 'proteingym-
                        dms/MET_HUMAN_Estevam_2023', 'proteingym-
                        dms/MK01_HUMAN_Brenan_2016', 'proteingym-
                        dms/MLAC_ECOLI_MacRae_2023', 'proteingym-
                        dms/MSH2_HUMAN_Jia_2020', 'proteingym-
                        dms/MTH3_HAEAE_RockahShmuel_2015', 'proteingym-
                        dms/MTHR_HUMAN_Weile_2021', 'proteingym-
                        dms/NCAP_I34A1_Doud_2015', 'proteingym-
                        dms/NPC1_HUMAN_Erwood_2022_HEK293T', 'proteingym-
                        dms/NPC1_HUMAN_Erwood_2022_RPE1', 'proteingym-
                        dms/NRAM_I33A0_Jiang_2016', 'proteingym-
                        dms/NUD15_HUMAN_Suiter_2020', 'proteingym-
                        dms/OPSD_HUMAN_Wan_2019', 'proteingym-
                        dms/OTC_HUMAN_Lo_2023', 'proteingym-
                        dms/OTU7A_HUMAN_Tsuboyama_2023_2L2D', 'proteingym-
                        dms/OXDA_RHOTO_Vanella_2023_activity', 'proteingym-
                        dms/OXDA_RHOTO_Vanella_2023_expression', 'proteingym-
                        dms/P53_HUMAN_Giacomelli_2018_Null_Etoposide',
                        'proteingym-
                        dms/P53_HUMAN_Giacomelli_2018_Null_Nutlin',
                        'proteingym-dms/P53_HUMAN_Giacomelli_2018_WT_Nutlin',
                        'proteingym-dms/P53_HUMAN_Kotler_2018', 'proteingym-
                        dms/P84126_THETH_Chan_2017', 'proteingym-
                        dms/PAI1_HUMAN_Huttinger_2021', 'proteingym-
                        dms/PA_I34A1_Wu_2015', 'proteingym-
                        dms/PKN1_HUMAN_Tsuboyama_2023_1URF', 'proteingym-
                        dms/POLG_CXB3N_Mattenberger_2021', 'proteingym-
                        dms/POLG_DEN26_Suphatrakul_2023', 'proteingym-
                        dms/POLG_HCVJF_Qi_2014', 'proteingym-
                        dms/PPARG_HUMAN_Majithia_2016', 'proteingym-
                        dms/PPM1D_HUMAN_Miller_2022', 'proteingym-
                        dms/PRKN_HUMAN_Clausen_2023', 'proteingym-
                        dms/PTEN_HUMAN_Matreyek_2021', 'proteingym-
                        dms/PTEN_HUMAN_Mighell_2018', 'proteingym-
                        dms/Q2N0S5_9HIV1_Haddox_2018', 'proteingym-
                        dms/Q53Z42_HUMAN_McShan_2019_binding-TAPBPR',
                        'proteingym-dms/Q53Z42_HUMAN_McShan_2019_expression',
                        'proteingym-dms/Q59976_STRSQ_Romero_2015',
                        'proteingym-dms/Q837P4_ENTFA_Meier_2023', 'proteingym-
                        dms/Q837P5_ENTFA_Meier_2023', 'proteingym-
                        dms/R1AB_SARS2_Flynn_2022', 'proteingym-
                        dms/RAF1_HUMAN_Zinkus-Boltz_2019', 'proteingym-
                        dms/RASH_HUMAN_Bandaru_2017', 'proteingym-
                        dms/RDRP_I33A0_Li_2023', 'proteingym-
                        dms/REV_HV1H2_Fernandes_2016', 'proteingym-
                        dms/RL40A_YEAST_Mavor_2016', 'proteingym-
                        dms/RL40A_YEAST_Roscoe_2013', 'proteingym-
                        dms/RL40A_YEAST_Roscoe_2014', 'proteingym-
                        dms/RNC_ECOLI_Weeks_2023', 'proteingym-
                        dms/RPC1_LAMBD_Li_2019_high-expression', 'proteingym-
                        dms/RPC1_LAMBD_Li_2019_low-expression', 'proteingym-
                        dms/RS15_GEOSE_Tsuboyama_2023_1A32', 'proteingym-
                        dms/S22A1_HUMAN_Yee_2023_abundance', 'proteingym-
                        dms/S22A1_HUMAN_Yee_2023_activity', 'proteingym-
                        dms/SBI_STAAM_Tsuboyama_2023_2JVG', 'proteingym-
                        dms/SC6A4_HUMAN_Young_2021', 'proteingym-
                        dms/SCIN_STAAR_Tsuboyama_2023_2QFF', 'proteingym-
                        dms/SCN5A_HUMAN_Glazer_2019', 'proteingym-
                        dms/SERC_HUMAN_Xie_2023', 'proteingym-
                        dms/SHOC2_HUMAN_Kwon_2022', 'proteingym-
                        dms/SOX30_HUMAN_Tsuboyama_2023_7JJK', 'proteingym-
                        dms/SPIKE_SARS2_Starr_2020_binding', 'proteingym-
                        dms/SPIKE_SARS2_Starr_2020_expression', 'proteingym-
                        dms/SQSTM_MOUSE_Tsuboyama_2023_2RRU', 'proteingym-
                        dms/SRC_HUMAN_Ahler_2019', 'proteingym-
                        dms/SRC_HUMAN_Chakraborty_2023_binding-DAS_25uM',
                        'proteingym-dms/SRC_HUMAN_Nguyen_2022', 'proteingym-
                        dms/SUMO1_HUMAN_Weile_2017', 'proteingym-
                        dms/SYUA_HUMAN_Newberry_2020', 'proteingym-
                        dms/TADBP_HUMAN_Bolognesi_2019', 'proteingym-
                        dms/TAT_HV1BR_Fernandes_2016', 'proteingym-
                        dms/TPK1_HUMAN_Weile_2017', 'proteingym-
                        dms/TPMT_HUMAN_Matreyek_2018', 'proteingym-
                        dms/TPOR_HUMAN_Bridgford_2020', 'proteingym-
                        dms/TRPC_SACS2_Chan_2017', 'proteingym-
                        dms/TRPC_THEMA_Chan_2017', 'proteingym-
                        dms/UBC9_HUMAN_Weile_2017', 'proteingym-
                        dms/UBE4B_MOUSE_Starita_2013', 'proteingym-
                        dms/VG08_BPP22_Tsuboyama_2023_2GP8', 'proteingym-
                        dms/VKOR1_HUMAN_Chiasson_2020_abundance', 'proteingym-
                        dms/VKOR1_HUMAN_Chiasson_2020_activity', 'proteingym-
                        dms/VRPI_BPT7_Tsuboyama_2023_2WNM']
  --all                 Download the leaderboard for all tasks.
  --force               Force download the leaderboard, even if it already
                        exists.

```

#### `biomlbench dev prepare-human-baselines`

```
usage: cli.py dev prepare-human-baselines [-h] [-t TASK_ID] [--all] [--force]

options:
  -h, --help            show this help message and exit
  -t TASK_ID, --task-id TASK_ID
                        Name of the task to prepare human baselines for. Valid
                        options: ['manual/caco2-wang',
                        'openproblems/batch_integration', 'polarishub/adaptyv-
                        bio-egfr-binders-binary-cls-v0', 'polarishub/biogen-
                        adme-fang-hclint-reg-v1', 'polarishub/biogen-adme-
                        fang-hppb-reg-v1', 'polarishub/biogen-adme-fang-perm-
                        reg-v1', 'polarishub/biogen-adme-fang-rclint-reg-v1',
                        'polarishub/biogen-adme-fang-reg-v1',
                        'polarishub/biogen-adme-fang-rppb-reg-v1',
                        'polarishub/biogen-adme-fang-solu-reg-v1',
                        'polarishub/graphium-l1000-mcf7-v1',
                        'polarishub/graphium-l1000-vcap-v1',
                        'polarishub/graphium-pcba-1328-1564k-v1',
                        'polarishub/graphium-qm9-v1', 'polarishub/graphium-
                        tox21-v1', 'polarishub/graphium-zinc12k-v1',
                        'polarishub/mlls-bend-zeroshot-variant-effects-
                        disease', 'polarishub/mlls-bend-zeroshot-variant-
                        effects-expression', 'polarishub/molecularml-
                        moleculeace-chembl1862-ki', 'polarishub/molecularml-
                        moleculeace-chembl1871-ki', 'polarishub/molecularml-
                        moleculeace-chembl2034-ki', 'polarishub/molecularml-
                        moleculeace-chembl204-ki', 'polarishub/molecularml-
                        moleculeace-chembl2047-ec50', 'polarishub/molecularml-
                        moleculeace-chembl214-ki', 'polarishub/molecularml-
                        moleculeace-chembl2147-ki', 'polarishub/molecularml-
                        moleculeace-chembl218-ec50', 'polarishub/molecularml-
                        moleculeace-chembl219-ki', 'polarishub/molecularml-
                        moleculeace-chembl228-ki', 'polarishub/molecularml-
                        moleculeace-chembl231-ki', 'polarishub/molecularml-
                        moleculeace-chembl233-ki', 'polarishub/molecularml-
                        moleculeace-chembl234-ki', 'polarishub/molecularml-
                        moleculeace-chembl235-ec50', 'polarishub/molecularml-
                        moleculeace-chembl236-ki', 'polarishub/molecularml-
                        moleculeace-chembl237-ec50', 'polarishub/molecularml-
                        moleculeace-chembl237-ki', 'polarishub/molecularml-
                        moleculeace-chembl238-ki', 'polarishub/molecularml-
                        moleculeace-chembl239-ec50', 'polarishub/molecularml-
                        moleculeace-chembl244-ki', 'polarishub/molecularml-
                        moleculeace-chembl262-ki', 'polarishub/molecularml-
                        moleculeace-chembl264-ki', 'polarishub/molecularml-
                        moleculeace-chembl2835-ki', 'polarishub/molecularml-
                        moleculeace-chembl287-ki', 'polarishub/molecularml-
                        moleculeace-chembl2971-ki', 'polarishub/molecularml-
                        moleculeace-chembl3979-ec50', 'polarishub/molecularml-
                        moleculeace-chembl4005-ki', 'polarishub/molecularml-
                        moleculeace-chembl4203-ki', 'polarishub/molecularml-
                        moleculeace-chembl4616-ec50', 'polarishub/molecularml-
                        moleculeace-chembl4792-ki', 'polarishub/novartis-adme-
                        novartis-cyp3a4-cls', 'polarishub/novartis-adme-
                        novartis-cyp3a4-reg', 'polarishub/polaris-adme-fang-
                        hclint-1', 'polarishub/polaris-adme-fang-hppb-1',
                        'polarishub/polaris-adme-fang-perm-1',
                        'polarishub/polaris-adme-fang-r-1',
                        'polarishub/polaris-adme-fang-rclint-1',
                        'polarishub/polaris-adme-fang-rppb-1',
                        'polarishub/polaris-adme-fang-solu-1',
                        'polarishub/polaris-hello-world-benchmark',
                        'polarishub/polaris-molprop-250k-leadlike-r-1',
                        'polarishub/polaris-molprop-250k-r-1',
                        'polarishub/polaris-molprop-250k-reg-v2',
                        'polarishub/polaris-molprop-leadlike-250k-reg-v2',
                        'polarishub/polaris-pkis1-egfr-wt-mut-c-1',
                        'polarishub/polaris-pkis1-egfr-wt-mut-r-1',
                        'polarishub/polaris-pkis1-kit-wt-mut-c-1',
                        'polarishub/polaris-pkis1-kit-wt-mut-r-1',
                        'polarishub/polaris-pkis1-ret-wt-mut-c-1',
                        'polarishub/polaris-pkis1-ret-wt-mut-r-1',
                        'polarishub/polaris-pkis2-egfr-wt-c-1',
                        'polarishub/polaris-pkis2-egfr-wt-r-1',
                        'polarishub/polaris-pkis2-egfr-wt-reg-v2',
                        'polarishub/polaris-pkis2-kit-wt-c-1',
                        'polarishub/polaris-pkis2-kit-wt-cls-v2',
                        'polarishub/polaris-pkis2-kit-wt-r-1',
                        'polarishub/polaris-pkis2-kit-wt-reg-v2',
                        'polarishub/polaris-pkis2-lok-slk-c-1',
                        'polarishub/polaris-pkis2-lok-slk-cls-v2',
                        'polarishub/polaris-pkis2-lok-slk-r-1',
                        'polarishub/polaris-pkis2-ret-wt-c-1',
                        'polarishub/polaris-pkis2-ret-wt-cls-v2',
                        'polarishub/polaris-pkis2-ret-wt-r-1',
                        'polarishub/polaris-pkis2-ret-wt-reg-v2',
                        'polarishub/polaris-posebusters-v1',
                        'polarishub/recursion-rxrx-compound-gene-activity-
                        benchmark', 'polarishub/tdcommons-ames',
                        'polarishub/tdcommons-bbb-martins',
                        'polarishub/tdcommons-bioavailability-ma',
                        'polarishub/tdcommons-caco2-wang',
                        'polarishub/tdcommons-clearance-hepatocyte-az',
                        'polarishub/tdcommons-clearance-microsome-az',
                        'polarishub/tdcommons-cyp2c9-substrate-carbonmangels',
                        'polarishub/tdcommons-cyp2c9-veith',
                        'polarishub/tdcommons-cyp2d6-substrate-carbonmangels',
                        'polarishub/tdcommons-cyp2d6-veith',
                        'polarishub/tdcommons-cyp3a4-substrate-carbonmangels',
                        'polarishub/tdcommons-cyp3a4-veith',
                        'polarishub/tdcommons-dili', 'polarishub/tdcommons-
                        half-life-obach', 'polarishub/tdcommons-herg',
                        'polarishub/tdcommons-hia-hou', 'polarishub/tdcommons-
                        ld50-zhu', 'polarishub/tdcommons-lipophilicity-
                        astrazeneca', 'polarishub/tdcommons-pgp-broccatelli',
                        'polarishub/tdcommons-ppbr-az', 'polarishub/tdcommons-
                        solubility-aqsoldb', 'polarishub/tdcommons-vdss-
                        lombardo', 'polarishub/vishrut64-cho-dna-expression-
                        prediction-dataset-task', 'polarishub/vishrut64-rna-
                        expression-prediction-dataset-task', 'proteingym-
                        dms/A0A140D2T1_ZIKV_Sourisseau_2019', 'proteingym-
                        dms/A0A192B1T2_9HIV1_Haddox_2018', 'proteingym-
                        dms/A0A1I9GEU1_NEIME_Kennouche_2019', 'proteingym-
                        dms/A0A247D711_LISMN_Stadelmann_2021', 'proteingym-
                        dms/A0A2Z5U3Z0_9INFA_Doud_2016', 'proteingym-
                        dms/A0A2Z5U3Z0_9INFA_Wu_2014', 'proteingym-
                        dms/A4D664_9INFA_Soh_2019', 'proteingym-
                        dms/A4GRB6_PSEAI_Chen_2020', 'proteingym-
                        dms/AACC1_PSEAI_Dandage_2018', 'proteingym-
                        dms/ACE2_HUMAN_Chan_2020', 'proteingym-
                        dms/ADRB2_HUMAN_Jones_2020', 'proteingym-
                        dms/AICDA_HUMAN_Gajula_2014_3cycles', 'proteingym-
                        dms/AMIE_PSEAE_Wrenbeck_2017', 'proteingym-
                        dms/ANCSZ_Hobbs_2022', 'proteingym-
                        dms/ARGR_ECOLI_Tsuboyama_2023_1AOY', 'proteingym-
                        dms/B2L11_HUMAN_Dutta_2010_binding-Mcl-1',
                        'proteingym-dms/BLAT_ECOLX_Deng_2012', 'proteingym-
                        dms/BLAT_ECOLX_Firnberg_2014', 'proteingym-
                        dms/BLAT_ECOLX_Jacquier_2013', 'proteingym-
                        dms/BLAT_ECOLX_Stiffler_2015', 'proteingym-
                        dms/BRCA1_HUMAN_Findlay_2018', 'proteingym-
                        dms/BRCA2_HUMAN_Erwood_2022_HEK293T', 'proteingym-
                        dms/C6KNH7_9INFA_Lee_2018', 'proteingym-
                        dms/CALM1_HUMAN_Weile_2017', 'proteingym-
                        dms/CAR11_HUMAN_Meitlis_2020_gof', 'proteingym-
                        dms/CAR11_HUMAN_Meitlis_2020_lof', 'proteingym-
                        dms/CAS9_STRP1_Spencer_2017_positive', 'proteingym-
                        dms/CASP3_HUMAN_Roychowdhury_2020', 'proteingym-
                        dms/CASP7_HUMAN_Roychowdhury_2020', 'proteingym-
                        dms/CBS_HUMAN_Sun_2020', 'proteingym-
                        dms/CCDB_ECOLI_Adkar_2012', 'proteingym-
                        dms/CCDB_ECOLI_Tripathi_2016', 'proteingym-
                        dms/CCR5_HUMAN_Gill_2023', 'proteingym-
                        dms/CD19_HUMAN_Klesmith_2019_FMC_singles',
                        'proteingym-dms/CP2C9_HUMAN_Amorosi_2021_abundance',
                        'proteingym-dms/CP2C9_HUMAN_Amorosi_2021_activity',
                        'proteingym-dms/DLG4_RAT_McLaughlin_2012',
                        'proteingym-dms/DN7A_SACS2_Tsuboyama_2023_1JIC',
                        'proteingym-dms/DYR_ECOLI_Nguyen_2023', 'proteingym-
                        dms/DYR_ECOLI_Thompson_2019', 'proteingym-
                        dms/ENVZ_ECOLI_Ghose_2023', 'proteingym-
                        dms/ENV_HV1B9_DuenasDecamp_2016', 'proteingym-
                        dms/ENV_HV1BR_Haddox_2016', 'proteingym-
                        dms/ERBB2_HUMAN_Elazar_2016', 'proteingym-
                        dms/ESTA_BACSU_Nutschel_2020', 'proteingym-
                        dms/FKBP3_HUMAN_Tsuboyama_2023_2KFV', 'proteingym-
                        dms/GAL4_YEAST_Kitzman_2015', 'proteingym-
                        dms/GDIA_HUMAN_Silverstein_2021', 'proteingym-
                        dms/GLPA_HUMAN_Elazar_2016', 'proteingym-
                        dms/HCP_LAMBD_Tsuboyama_2023_2L6Q', 'proteingym-
                        dms/HEM3_HUMAN_Loggerenberg_2023', 'proteingym-
                        dms/HMDH_HUMAN_Jiang_2019', 'proteingym-
                        dms/HSP82_YEAST_Cote-Hammarlof_2020_growth-H2O2',
                        'proteingym-dms/HSP82_YEAST_Flynn_2019', 'proteingym-
                        dms/HSP82_YEAST_Mishra_2016', 'proteingym-
                        dms/HXK4_HUMAN_Gersing_2022_activity', 'proteingym-
                        dms/HXK4_HUMAN_Gersing_2023_abundance', 'proteingym-
                        dms/I6TAH8_I68A0_Doud_2015', 'proteingym-
                        dms/IF1_ECOLI_Kelsic_2016', 'proteingym-
                        dms/ILF3_HUMAN_Tsuboyama_2023_2L33', 'proteingym-
                        dms/KCNE1_HUMAN_Muhammad_2023_expression',
                        'proteingym-dms/KCNE1_HUMAN_Muhammad_2023_function',
                        'proteingym-dms/KCNH2_HUMAN_Kozek_2020', 'proteingym-
                        dms/KCNJ2_MOUSE_Coyote-Maestas_2022_function',
                        'proteingym-dms/KCNJ2_MOUSE_Coyote-
                        Maestas_2022_surface', 'proteingym-
                        dms/KKA2_KLEPN_Melnikov_2014', 'proteingym-
                        dms/LGK_LIPST_Klesmith_2015', 'proteingym-
                        dms/LYAM1_HUMAN_Elazar_2016', 'proteingym-
                        dms/MET_HUMAN_Estevam_2023', 'proteingym-
                        dms/MK01_HUMAN_Brenan_2016', 'proteingym-
                        dms/MLAC_ECOLI_MacRae_2023', 'proteingym-
                        dms/MSH2_HUMAN_Jia_2020', 'proteingym-
                        dms/MTH3_HAEAE_RockahShmuel_2015', 'proteingym-
                        dms/MTHR_HUMAN_Weile_2021', 'proteingym-
                        dms/NCAP_I34A1_Doud_2015', 'proteingym-
                        dms/NPC1_HUMAN_Erwood_2022_HEK293T', 'proteingym-
                        dms/NPC1_HUMAN_Erwood_2022_RPE1', 'proteingym-
                        dms/NRAM_I33A0_Jiang_2016', 'proteingym-
                        dms/NUD15_HUMAN_Suiter_2020', 'proteingym-
                        dms/OPSD_HUMAN_Wan_2019', 'proteingym-
                        dms/OTC_HUMAN_Lo_2023', 'proteingym-
                        dms/OTU7A_HUMAN_Tsuboyama_2023_2L2D', 'proteingym-
                        dms/OXDA_RHOTO_Vanella_2023_activity', 'proteingym-
                        dms/OXDA_RHOTO_Vanella_2023_expression', 'proteingym-
                        dms/P53_HUMAN_Giacomelli_2018_Null_Etoposide',
                        'proteingym-
                        dms/P53_HUMAN_Giacomelli_2018_Null_Nutlin',
                        'proteingym-dms/P53_HUMAN_Giacomelli_2018_WT_Nutlin',
                        'proteingym-dms/P53_HUMAN_Kotler_2018', 'proteingym-
                        dms/P84126_THETH_Chan_2017', 'proteingym-
                        dms/PAI1_HUMAN_Huttinger_2021', 'proteingym-
                        dms/PA_I34A1_Wu_2015', 'proteingym-
                        dms/PKN1_HUMAN_Tsuboyama_2023_1URF', 'proteingym-
                        dms/POLG_CXB3N_Mattenberger_2021', 'proteingym-
                        dms/POLG_DEN26_Suphatrakul_2023', 'proteingym-
                        dms/POLG_HCVJF_Qi_2014', 'proteingym-
                        dms/PPARG_HUMAN_Majithia_2016', 'proteingym-
                        dms/PPM1D_HUMAN_Miller_2022', 'proteingym-
                        dms/PRKN_HUMAN_Clausen_2023', 'proteingym-
                        dms/PTEN_HUMAN_Matreyek_2021', 'proteingym-
                        dms/PTEN_HUMAN_Mighell_2018', 'proteingym-
                        dms/Q2N0S5_9HIV1_Haddox_2018', 'proteingym-
                        dms/Q53Z42_HUMAN_McShan_2019_binding-TAPBPR',
                        'proteingym-dms/Q53Z42_HUMAN_McShan_2019_expression',
                        'proteingym-dms/Q59976_STRSQ_Romero_2015',
                        'proteingym-dms/Q837P4_ENTFA_Meier_2023', 'proteingym-
                        dms/Q837P5_ENTFA_Meier_2023', 'proteingym-
                        dms/R1AB_SARS2_Flynn_2022', 'proteingym-
                        dms/RAF1_HUMAN_Zinkus-Boltz_2019', 'proteingym-
                        dms/RASH_HUMAN_Bandaru_2017', 'proteingym-
                        dms/RDRP_I33A0_Li_2023', 'proteingym-
                        dms/REV_HV1H2_Fernandes_2016', 'proteingym-
                        dms/RL40A_YEAST_Mavor_2016', 'proteingym-
                        dms/RL40A_YEAST_Roscoe_2013', 'proteingym-
                        dms/RL40A_YEAST_Roscoe_2014', 'proteingym-
                        dms/RNC_ECOLI_Weeks_2023', 'proteingym-
                        dms/RPC1_LAMBD_Li_2019_high-expression', 'proteingym-
                        dms/RPC1_LAMBD_Li_2019_low-expression', 'proteingym-
                        dms/RS15_GEOSE_Tsuboyama_2023_1A32', 'proteingym-
                        dms/S22A1_HUMAN_Yee_2023_abundance', 'proteingym-
                        dms/S22A1_HUMAN_Yee_2023_activity', 'proteingym-
                        dms/SBI_STAAM_Tsuboyama_2023_2JVG', 'proteingym-
                        dms/SC6A4_HUMAN_Young_2021', 'proteingym-
                        dms/SCIN_STAAR_Tsuboyama_2023_2QFF', 'proteingym-
                        dms/SCN5A_HUMAN_Glazer_2019', 'proteingym-
                        dms/SERC_HUMAN_Xie_2023', 'proteingym-
                        dms/SHOC2_HUMAN_Kwon_2022', 'proteingym-
                        dms/SOX30_HUMAN_Tsuboyama_2023_7JJK', 'proteingym-
                        dms/SPIKE_SARS2_Starr_2020_binding', 'proteingym-
                        dms/SPIKE_SARS2_Starr_2020_expression', 'proteingym-
                        dms/SQSTM_MOUSE_Tsuboyama_2023_2RRU', 'proteingym-
                        dms/SRC_HUMAN_Ahler_2019', 'proteingym-
                        dms/SRC_HUMAN_Chakraborty_2023_binding-DAS_25uM',
                        'proteingym-dms/SRC_HUMAN_Nguyen_2022', 'proteingym-
                        dms/SUMO1_HUMAN_Weile_2017', 'proteingym-
                        dms/SYUA_HUMAN_Newberry_2020', 'proteingym-
                        dms/TADBP_HUMAN_Bolognesi_2019', 'proteingym-
                        dms/TAT_HV1BR_Fernandes_2016', 'proteingym-
                        dms/TPK1_HUMAN_Weile_2017', 'proteingym-
                        dms/TPMT_HUMAN_Matreyek_2018', 'proteingym-
                        dms/TPOR_HUMAN_Bridgford_2020', 'proteingym-
                        dms/TRPC_SACS2_Chan_2017', 'proteingym-
                        dms/TRPC_THEMA_Chan_2017', 'proteingym-
                        dms/UBC9_HUMAN_Weile_2017', 'proteingym-
                        dms/UBE4B_MOUSE_Starita_2013', 'proteingym-
                        dms/VG08_BPP22_Tsuboyama_2023_2GP8', 'proteingym-
                        dms/VKOR1_HUMAN_Chiasson_2020_abundance', 'proteingym-
                        dms/VKOR1_HUMAN_Chiasson_2020_activity', 'proteingym-
                        dms/VRPI_BPT7_Tsuboyama_2023_2WNM']
  --all                 Prepare human baselines for all tasks.
  --force               Force re-extraction of human baselines, even if they
                        already exist.

```

## Usage Examples

### Task Preparation

```bash
# Prepare a specific task
biomlbench prepare -t polarishub/tdcommons-caco2-wang

# Prepare all tasks in a domain
biomlbench prepare --domain admet

# Prepare multiple tasks from a file
biomlbench prepare --list experiments/splits/caco2-wang.txt

# Prepare all low-difficulty tasks
biomlbench prepare --lite
```

### Agent Execution

```bash
# Run agent on single task
biomlbench run-agent --agent dummy --task-id polarishub/tdcommons-caco2-wang

# Run agent on multiple tasks with parallel workers
biomlbench run-agent \
    --agent aide \
    --task-list experiments/splits/caco2-wang.txt \
    --n-workers 4 \
    --n-seeds 3

# Run with custom container configuration
biomlbench run-agent \
    --agent aide \
    --task-id polarishub/tdcommons-caco2-wang \
    --container-config custom_config.json \
    --retain-container
```

### Evaluation and Grading

```bash
# Grade multiple task submissions
biomlbench grade \
    --submission runs/my-run-group/submission.jsonl \
    --output-dir results/

# Grade single task submission (CSV format)
biomlbench grade-sample submission.csv polarishub/tdcommons-caco2-wang

# Grade single task submission (H5AD format)
biomlbench grade-sample submission.h5ad openproblems/cell_cell_communication

# Run and grade baselines
biomlbench run-baseline polarishub/tdcommons-caco2-wang --baseline all
biomlbench grade \
    --submission baseline_submissions/submission.jsonl \
    --output-dir baseline_results/
```

### Development Commands

```bash
# Download leaderboards for all tasks
biomlbench dev download-leaderboard --all

# Prepare human baselines for a specific task
biomlbench dev prepare-human-baselines -t polarishub/tdcommons-caco2-wang --force
```

## Environment Variables

BioML-bench respects these environment variables:

### Agent Configuration  
- **`OPENAI_API_KEY`** - API key for AIDE agent
- **`I_ACCEPT_RUNNING_PRIVILEGED_CONTAINERS`** - Set to "true" to allow privileged containers

### Kaggle Integration
Kaggle authentication uses the standard configuration file at `~/.kaggle/kaggle.json`. See the [Kaggle API documentation](https://github.com/Kaggle/kaggle-api) for setup instructions.

## Common Patterns

### Full Workflow Example

```bash
# 1. Prepare tasks
biomlbench prepare --lite

# 2. Run agent
biomlbench run-agent --agent dummy --task-list experiments/splits/polarishub-tdcommons-caco2-wang.txt

# 3. Grade results (submission.jsonl is auto-generated)
biomlbench grade \
    --submission runs/latest-run-group/submission.jsonl \
    --output-dir results/
```

### Debugging Agent Issues

```bash
# Run with container retention for debugging
biomlbench run-agent \
    --agent my-agent \
    --task-id polarishub/tdcommons-caco2-wang \
    --retain-container

# Check logs in the run directory
ls runs/latest-run-group/*/logs/
```
