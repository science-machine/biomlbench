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
    grade-sample        Grade a single sample (task) in the eval
    dev                 Developer tools for extending BioML-bench.
    run-baseline        Run a baseline agent on a biomedical task
    run-agent           Run an AI agent on biomedical tasks

options:
  -h, --help            show this help message and exit

```

## Commands

### `biomlbench prepare`

```
usage: cli.py prepare [-h] [-t TASK_ID] [-a] [--lite] [-l LIST]
                      [--domain DOMAIN] [--task-type TASK_TYPE] [--keep-raw]
                      [--data-dir DATA_DIR] [--overwrite-checksums]
                      [--overwrite-leaderboard] [--skip-verification]

options:
  -h, --help            show this help message and exit
  -t TASK_ID, --task-id TASK_ID
                        ID of the task to prepare in 'folder/task' format.
                        Examples: manual/caco2-wang, polarishub/tdcommons-
                        admet
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
                        [--output-dir OUTPUT_DIR] [--s3-bucket S3_BUCKET]
                        [--s3-prefix S3_PREFIX] [--s3-no-compress]
                        [--disable-s3] [--incremental-s3] [--cpu-only]
                        [--fast]

options:
  -h, --help            show this help message and exit
  --agent AGENT         Agent ID to run (e.g., dummy, aide, aide/dev)
  --task-id TASK_ID     Single task ID to run. Valid options:
                        ['kaggle/histopathologic-cancer-detection',
                        'kaggle/hms-harmful-brain-activity-classification',
                        'kaggle/osic-pulmonary-fibrosis-progression',
                        'kaggle/ranzcr-clip-catheter-line-classification',
                        'kaggle/rsna-miccai-brain-tumor-radiogenomic-
                        classification', 'kaggle/stanford-covid-vaccine',
                        'kaggle/uw-madison-gi-tract-image-segmentation',
                        'kaggle/ventilator-pressure-prediction',
                        'kaggle_big/bms-molecular-translation',
                        'kaggle_big/rsna-2022-cervical-spine-fracture-
                        detection', 'kaggle_big/rsna-breast-cancer-detection',
                        'kaggle_big/siim-covid19-detection', 'kaggle_big/siim-
                        isic-melanoma-classification', 'kaggle_big/vinbigdata-
                        chest-xray-abnormalities-detection', 'manual/open-
                        problems-cell-cell-communication-ligand-target',
                        'manual/open-problems-label-projection', 'manual/open-
                        problems-predict-modality', 'manual/open-problems-
                        single-cell-perturbations', 'manual/open-problems-
                        spatially-variable-genes', 'polarishub/polaris-adme-
                        fang-hclint-1', 'polarishub/polaris-adme-fang-hppb-1',
                        'polarishub/polaris-adme-fang-perm-1',
                        'polarishub/polaris-adme-fang-r-1',
                        'polarishub/polaris-adme-fang-rclint-1',
                        'polarishub/polaris-adme-fang-rppb-1',
                        'polarishub/polaris-adme-fang-solu-1',
                        'polarishub/polaris-pkis1-egfr-wt-mut-c-1',
                        'polarishub/polaris-pkis1-egfr-wt-mut-r-1',
                        'polarishub/polaris-pkis1-kit-wt-mut-c-1',
                        'polarishub/polaris-pkis1-kit-wt-mut-r-1',
                        'polarishub/polaris-pkis1-ret-wt-mut-c-1',
                        'polarishub/polaris-pkis1-ret-wt-mut-r-1',
                        'polarishub/polaris-pkis2-egfr-wt-c-1',
                        'polarishub/polaris-pkis2-egfr-wt-r-1',
                        'polarishub/polaris-pkis2-kit-wt-c-1',
                        'polarishub/polaris-pkis2-kit-wt-r-1',
                        'polarishub/polaris-pkis2-lok-slk-c-1',
                        'polarishub/polaris-pkis2-lok-slk-r-1',
                        'polarishub/polaris-pkis2-ret-wt-c-1',
                        'polarishub/polaris-pkis2-ret-wt-r-1',
                        'polarishub/tdcommons-bbb-martins',
                        'polarishub/tdcommons-bioavailability-ma',
                        'polarishub/tdcommons-caco2-wang',
                        'polarishub/tdcommons-clearance-hepatocyte-az',
                        'polarishub/tdcommons-cyp2d6-substrate-carbonmangels',
                        'polarishub/tdcommons-herg', 'polarishub/tdcommons-
                        hia-hou', 'polarishub/tdcommons-lipophilicity-
                        astrazeneca', 'polarishub/tdcommons-pgp-broccatelli',
                        'polarishub-extended/adaptyv-bio-egfr-binders-binary-
                        cls-v0', 'polarishub-extended/biogen-adme-fang-hclint-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-hppb-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-perm-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-rclint-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-rppb-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-solu-
                        reg-v1', 'polarishub-extended/graphium-l1000-mcf7-v1',
                        'polarishub-extended/graphium-l1000-vcap-v1',
                        'polarishub-extended/graphium-pcba-1328-1564k-v1',
                        'polarishub-extended/graphium-qm9-v1', 'polarishub-
                        extended/graphium-tox21-v1', 'polarishub-
                        extended/graphium-zinc12k-v1', 'polarishub-
                        extended/novartis-adme-novartis-cyp3a4-cls',
                        'polarishub-extended/novartis-adme-novartis-
                        cyp3a4-reg', 'polarishub-extended/polaris-adme-fang-
                        hclint-1', 'polarishub-extended/polaris-adme-fang-
                        hppb-1', 'polarishub-extended/polaris-adme-fang-
                        perm-1', 'polarishub-extended/polaris-adme-fang-r-1',
                        'polarishub-extended/polaris-adme-fang-rclint-1',
                        'polarishub-extended/polaris-adme-fang-rppb-1',
                        'polarishub-extended/polaris-adme-fang-solu-1',
                        'polarishub-extended/polaris-hello-world-benchmark',
                        'polarishub-extended/polaris-
                        molprop-250k-leadlike-r-1', 'polarishub-
                        extended/polaris-molprop-250k-r-1', 'polarishub-
                        extended/polaris-molprop-250k-reg-v2', 'polarishub-
                        extended/polaris-molprop-leadlike-250k-reg-v2',
                        'polarishub-extended/polaris-pkis1-egfr-wt-mut-c-1',
                        'polarishub-extended/polaris-pkis1-egfr-wt-mut-r-1',
                        'polarishub-extended/polaris-pkis1-kit-wt-mut-c-1',
                        'polarishub-extended/polaris-pkis1-kit-wt-mut-r-1',
                        'polarishub-extended/polaris-pkis1-ret-wt-mut-c-1',
                        'polarishub-extended/polaris-pkis1-ret-wt-mut-r-1',
                        'polarishub-extended/polaris-pkis2-egfr-wt-c-1',
                        'polarishub-extended/polaris-pkis2-egfr-wt-r-1',
                        'polarishub-extended/polaris-pkis2-egfr-wt-reg-v2',
                        'polarishub-extended/polaris-pkis2-kit-wt-c-1',
                        'polarishub-extended/polaris-pkis2-kit-wt-cls-v2',
                        'polarishub-extended/polaris-pkis2-kit-wt-r-1',
                        'polarishub-extended/polaris-pkis2-kit-wt-reg-v2',
                        'polarishub-extended/polaris-pkis2-lok-slk-c-1',
                        'polarishub-extended/polaris-pkis2-lok-slk-cls-v2',
                        'polarishub-extended/polaris-pkis2-lok-slk-r-1',
                        'polarishub-extended/polaris-pkis2-ret-wt-c-1',
                        'polarishub-extended/polaris-pkis2-ret-wt-cls-v2',
                        'polarishub-extended/polaris-pkis2-ret-wt-r-1',
                        'polarishub-extended/polaris-pkis2-ret-wt-reg-v2',
                        'polarishub-extended/polaris-posebusters-v1',
                        'polarishub-extended/recursion-rxrx-compound-gene-
                        activity-benchmark', 'polarishub-extended/tdcommons-
                        ames', 'polarishub-extended/tdcommons-bbb-martins',
                        'polarishub-extended/tdcommons-bioavailability-ma',
                        'polarishub-extended/tdcommons-caco2-wang',
                        'polarishub-extended/tdcommons-clearance-hepatocyte-
                        az', 'polarishub-extended/tdcommons-clearance-
                        microsome-az', 'polarishub-extended/tdcommons-
                        cyp2c9-substrate-carbonmangels', 'polarishub-
                        extended/tdcommons-cyp2c9-veith', 'polarishub-
                        extended/tdcommons-cyp2d6-substrate-carbonmangels',
                        'polarishub-extended/tdcommons-cyp2d6-veith',
                        'polarishub-extended/tdcommons-cyp3a4-substrate-
                        carbonmangels', 'polarishub-extended/tdcommons-
                        cyp3a4-veith', 'polarishub-extended/tdcommons-dili',
                        'polarishub-extended/tdcommons-half-life-obach',
                        'polarishub-extended/tdcommons-herg', 'polarishub-
                        extended/tdcommons-hia-hou', 'polarishub-
                        extended/tdcommons-ld50-zhu', 'polarishub-
                        extended/tdcommons-lipophilicity-astrazeneca',
                        'polarishub-extended/tdcommons-pgp-broccatelli',
                        'polarishub-extended/tdcommons-ppbr-az', 'polarishub-
                        extended/tdcommons-solubility-aqsoldb', 'polarishub-
                        extended/tdcommons-vdss-lombardo', 'polarishub-
                        extended/vishrut64-cho-dna-expression-prediction-
                        dataset-task', 'polarishub-extended/vishrut64-rna-
                        expression-prediction-dataset-task', 'proteingym-
                        dms/A0A140D2T1_ZIKV_Sourisseau_2019', 'proteingym-
                        dms/A0A192B1T2_9HIV1_Haddox_2018', 'proteingym-
                        dms/A0A1I9GEU1_NEIME_Kennouche_2019', 'proteingym-
                        dms/A0A247D711_LISMN_Stadelmann_2021', 'proteingym-
                        dms/A0A2Z5U3Z0_9INFA_Doud_2016', 'proteingym-
                        dms/A0A2Z5U3Z0_9INFA_Wu_2014', 'proteingym-
                        dms/A4D664_9INFA_Soh_2019', 'proteingym-
                        dms/A4GRB6_PSEAI_Chen_2020', 'proteingym-
                        dms/A4_HUMAN_Seuma_2022', 'proteingym-
                        dms/A4_HUMAN_Seuma_2022_indels', 'proteingym-
                        dms/AACC1_PSEAI_Dandage_2018', 'proteingym-
                        dms/ACE2_HUMAN_Chan_2020', 'proteingym-
                        dms/ADRB2_HUMAN_Jones_2020', 'proteingym-
                        dms/AICDA_HUMAN_Gajula_2014_3cycles', 'proteingym-
                        dms/AMFR_HUMAN_Tsuboyama_2023_4G3O', 'proteingym-
                        dms/AMFR_HUMAN_Tsuboyama_2023_4G3O_indels',
                        'proteingym-dms/AMIE_PSEAE_Wrenbeck_2017',
                        'proteingym-dms/ANCSZ_Hobbs_2022', 'proteingym-
                        dms/ARGR_ECOLI_Tsuboyama_2023_1AOY', 'proteingym-
                        dms/ARGR_ECOLI_Tsuboyama_2023_1AOY_indels',
                        'proteingym-dms/B1LPA6_ECOSM_Russ_2020_indels',
                        'proteingym-dms/B2L11_HUMAN_Dutta_2010_binding-Mcl-1',
                        'proteingym-dms/BBC1_YEAST_Tsuboyama_2023_1TG0',
                        'proteingym-
                        dms/BBC1_YEAST_Tsuboyama_2023_1TG0_indels',
                        'proteingym-dms/BCHB_CHLTE_Tsuboyama_2023_2KRU',
                        'proteingym-
                        dms/BCHB_CHLTE_Tsuboyama_2023_2KRU_indels',
                        'proteingym-dms/BLAT_ECOLX_Deng_2012', 'proteingym-
                        dms/BLAT_ECOLX_Firnberg_2014', 'proteingym-
                        dms/BLAT_ECOLX_Gonzalez_2019_indels', 'proteingym-
                        dms/BLAT_ECOLX_Jacquier_2013', 'proteingym-
                        dms/BLAT_ECOLX_Stiffler_2015', 'proteingym-
                        dms/BRCA1_HUMAN_Findlay_2018', 'proteingym-
                        dms/BRCA2_HUMAN_Erwood_2022_HEK293T', 'proteingym-
                        dms/C6KNH7_9INFA_Lee_2018', 'proteingym-
                        dms/CALM1_HUMAN_Weile_2017', 'proteingym-
                        dms/CAPSD_AAV2S_Sinai_2021', 'proteingym-
                        dms/CAPSD_AAV2S_Sinai_2021_designed_indels',
                        'proteingym-
                        dms/CAPSD_AAV2S_Sinai_2021_library_indels',
                        'proteingym-dms/CAR11_HUMAN_Meitlis_2020_gof',
                        'proteingym-dms/CAR11_HUMAN_Meitlis_2020_lof',
                        'proteingym-dms/CAS9_STRP1_Spencer_2017_positive',
                        'proteingym-dms/CASP3_HUMAN_Roychowdhury_2020',
                        'proteingym-dms/CASP7_HUMAN_Roychowdhury_2020',
                        'proteingym-dms/CATR_CHLRE_Tsuboyama_2023_2AMI',
                        'proteingym-
                        dms/CATR_CHLRE_Tsuboyama_2023_2AMI_indels',
                        'proteingym-dms/CBPA2_HUMAN_Tsuboyama_2023_1O6X',
                        'proteingym-
                        dms/CBPA2_HUMAN_Tsuboyama_2023_1O6X_indels',
                        'proteingym-dms/CBS_HUMAN_Sun_2020', 'proteingym-
                        dms/CBX4_HUMAN_Tsuboyama_2023_2K28', 'proteingym-
                        dms/CBX4_HUMAN_Tsuboyama_2023_2K28_indels',
                        'proteingym-dms/CCDB_ECOLI_Adkar_2012', 'proteingym-
                        dms/CCDB_ECOLI_Tripathi_2016', 'proteingym-
                        dms/CCR5_HUMAN_Gill_2023', 'proteingym-
                        dms/CD19_HUMAN_Klesmith_2019_FMC_singles',
                        'proteingym-dms/CP2C9_HUMAN_Amorosi_2021_abundance',
                        'proteingym-dms/CP2C9_HUMAN_Amorosi_2021_activity',
                        'proteingym-dms/CSN4_MOUSE_Tsuboyama_2023_1UFM',
                        'proteingym-
                        dms/CSN4_MOUSE_Tsuboyama_2023_1UFM_indels',
                        'proteingym-dms/CUE1_YEAST_Tsuboyama_2023_2MYX',
                        'proteingym-
                        dms/CUE1_YEAST_Tsuboyama_2023_2MYX_indels',
                        'proteingym-dms/D7PM05_CLYGR_Somermeyer_2022',
                        'proteingym-dms/DLG4_HUMAN_Faure_2021', 'proteingym-
                        dms/DLG4_RAT_McLaughlin_2012', 'proteingym-
                        dms/DN7A_SACS2_Tsuboyama_2023_1JIC', 'proteingym-
                        dms/DN7A_SACS2_Tsuboyama_2023_1JIC_indels',
                        'proteingym-dms/DNJA1_HUMAN_Tsuboyama_2023_2LO1',
                        'proteingym-
                        dms/DNJA1_HUMAN_Tsuboyama_2023_2LO1_indels',
                        'proteingym-dms/DOCK1_MOUSE_Tsuboyama_2023_2M0Y',
                        'proteingym-
                        dms/DOCK1_MOUSE_Tsuboyama_2023_2M0Y_indels',
                        'proteingym-dms/DYR_ECOLI_Nguyen_2023', 'proteingym-
                        dms/DYR_ECOLI_Thompson_2019', 'proteingym-
                        dms/ENVZ_ECOLI_Ghose_2023', 'proteingym-
                        dms/ENV_HV1B9_DuenasDecamp_2016', 'proteingym-
                        dms/ENV_HV1BR_Haddox_2016', 'proteingym-
                        dms/EPHB2_HUMAN_Tsuboyama_2023_1F0M', 'proteingym-
                        dms/EPHB2_HUMAN_Tsuboyama_2023_1F0M_indels',
                        'proteingym-dms/ERBB2_HUMAN_Elazar_2016', 'proteingym-
                        dms/ESTA_BACSU_Nutschel_2020', 'proteingym-
                        dms/F7YBW8_MESOW_Aakre_2015', 'proteingym-
                        dms/F7YBW8_MESOW_Ding_2023', 'proteingym-
                        dms/FECA_ECOLI_Tsuboyama_2023_2D1U', 'proteingym-
                        dms/FECA_ECOLI_Tsuboyama_2023_2D1U_indels',
                        'proteingym-dms/FKBP3_HUMAN_Tsuboyama_2023_2KFV',
                        'proteingym-dms/GAL4_YEAST_Kitzman_2015', 'proteingym-
                        dms/GCN4_YEAST_Staller_2018', 'proteingym-
                        dms/GDIA_HUMAN_Silverstein_2021', 'proteingym-
                        dms/GFP_AEQVI_Sarkisyan_2016', 'proteingym-
                        dms/GLPA_HUMAN_Elazar_2016', 'proteingym-
                        dms/GRB2_HUMAN_Faure_2021', 'proteingym-
                        dms/HCP_LAMBD_Tsuboyama_2023_2L6Q', 'proteingym-
                        dms/HCP_LAMBD_Tsuboyama_2023_2L6Q_indels',
                        'proteingym-dms/HECD1_HUMAN_Tsuboyama_2023_3DKM',
                        'proteingym-
                        dms/HECD1_HUMAN_Tsuboyama_2023_3DKM_indels',
                        'proteingym-dms/HEM3_HUMAN_Loggerenberg_2023',
                        'proteingym-dms/HIS7_YEAST_Pokusaeva_2019',
                        'proteingym-dms/HIS7_YEAST_Pokusaeva_2019_indels',
                        'proteingym-dms/HMDH_HUMAN_Jiang_2019', 'proteingym-
                        dms/HSP82_YEAST_Cote-Hammarlof_2020_growth-H2O2',
                        'proteingym-dms/HSP82_YEAST_Flynn_2019', 'proteingym-
                        dms/HSP82_YEAST_Mishra_2016', 'proteingym-
                        dms/HXK4_HUMAN_Gersing_2022_activity', 'proteingym-
                        dms/HXK4_HUMAN_Gersing_2023_abundance', 'proteingym-
                        dms/I6TAH8_I68A0_Doud_2015', 'proteingym-
                        dms/IF1_ECOLI_Kelsic_2016', 'proteingym-
                        dms/ILF3_HUMAN_Tsuboyama_2023_2L33', 'proteingym-
                        dms/ILF3_HUMAN_Tsuboyama_2023_2L33_indels',
                        'proteingym-dms/ISDH_STAAW_Tsuboyama_2023_2LHR',
                        'proteingym-dms/KCNE1_HUMAN_Muhammad_2023_expression',
                        'proteingym-dms/KCNE1_HUMAN_Muhammad_2023_function',
                        'proteingym-dms/KCNH2_HUMAN_Kozek_2020', 'proteingym-
                        dms/KCNJ2_MOUSE_Coyote-Maestas_2022_function',
                        'proteingym-dms/KCNJ2_MOUSE_Coyote-
                        Maestas_2022_surface', 'proteingym-
                        dms/KCNJ2_MOUSE_Macdonald_2022_indels', 'proteingym-
                        dms/KKA2_KLEPN_Melnikov_2014', 'proteingym-
                        dms/LGK_LIPST_Klesmith_2015', 'proteingym-
                        dms/LYAM1_HUMAN_Elazar_2016', 'proteingym-
                        dms/MAFG_MOUSE_Tsuboyama_2023_1K1V', 'proteingym-
                        dms/MAFG_MOUSE_Tsuboyama_2023_1K1V_indels',
                        'proteingym-dms/MBD11_ARATH_Tsuboyama_2023_6ACV',
                        'proteingym-
                        dms/MBD11_ARATH_Tsuboyama_2023_6ACV_indels',
                        'proteingym-dms/MET_HUMAN_Estevam_2023', 'proteingym-
                        dms/MK01_HUMAN_Brenan_2016', 'proteingym-
                        dms/MLAC_ECOLI_MacRae_2023', 'proteingym-
                        dms/MSH2_HUMAN_Jia_2020', 'proteingym-
                        dms/MTH3_HAEAE_RockahShmuel_2015', 'proteingym-
                        dms/MTHR_HUMAN_Weile_2021', 'proteingym-
                        dms/MYO3_YEAST_Tsuboyama_2023_2BTT', 'proteingym-
                        dms/MYO3_YEAST_Tsuboyama_2023_2BTT_indels',
                        'proteingym-dms/NCAP_I34A1_Doud_2015', 'proteingym-
                        dms/NKX31_HUMAN_Tsuboyama_2023_2L9R', 'proteingym-
                        dms/NKX31_HUMAN_Tsuboyama_2023_2L9R_indels',
                        'proteingym-dms/NPC1_HUMAN_Erwood_2022_HEK293T',
                        'proteingym-dms/NPC1_HUMAN_Erwood_2022_RPE1',
                        'proteingym-dms/NRAM_I33A0_Jiang_2016', 'proteingym-
                        dms/NUD15_HUMAN_Suiter_2020', 'proteingym-
                        dms/NUSA_ECOLI_Tsuboyama_2023_1WCL', 'proteingym-
                        dms/NUSA_ECOLI_Tsuboyama_2023_1WCL_indels',
                        'proteingym-dms/NUSG_MYCTU_Tsuboyama_2023_2MI6',
                        'proteingym-
                        dms/NUSG_MYCTU_Tsuboyama_2023_2MI6_indels',
                        'proteingym-dms/OBSCN_HUMAN_Tsuboyama_2023_1V1C',
                        'proteingym-
                        dms/OBSCN_HUMAN_Tsuboyama_2023_1V1C_indels',
                        'proteingym-dms/ODP2_GEOSE_Tsuboyama_2023_1W4G',
                        'proteingym-
                        dms/ODP2_GEOSE_Tsuboyama_2023_1W4G_indels',
                        'proteingym-dms/OPSD_HUMAN_Wan_2019', 'proteingym-
                        dms/OTC_HUMAN_Lo_2023', 'proteingym-
                        dms/OTU7A_HUMAN_Tsuboyama_2023_2L2D', 'proteingym-
                        dms/OTU7A_HUMAN_Tsuboyama_2023_2L2D_indels',
                        'proteingym-dms/OXDA_RHOTO_Vanella_2023_activity',
                        'proteingym-dms/OXDA_RHOTO_Vanella_2023_expression',
                        'proteingym-
                        dms/P53_HUMAN_Giacomelli_2018_Null_Etoposide',
                        'proteingym-
                        dms/P53_HUMAN_Giacomelli_2018_Null_Nutlin',
                        'proteingym-dms/P53_HUMAN_Giacomelli_2018_WT_Nutlin',
                        'proteingym-dms/P53_HUMAN_Kotler_2018', 'proteingym-
                        dms/P53_HUMAN_Kotler_2018_indels', 'proteingym-
                        dms/P84126_THETH_Chan_2017', 'proteingym-
                        dms/PABP_YEAST_Melamed_2013', 'proteingym-
                        dms/PAI1_HUMAN_Huttinger_2021', 'proteingym-
                        dms/PA_I34A1_Wu_2015', 'proteingym-
                        dms/PHOT_CHLRE_Chen_2023', 'proteingym-
                        dms/PIN1_HUMAN_Tsuboyama_2023_1I6C', 'proteingym-
                        dms/PIN1_HUMAN_Tsuboyama_2023_1I6C_indels',
                        'proteingym-dms/PITX2_HUMAN_Tsuboyama_2023_2L7M',
                        'proteingym-
                        dms/PITX2_HUMAN_Tsuboyama_2023_2L7M_indels',
                        'proteingym-dms/PKN1_HUMAN_Tsuboyama_2023_1URF',
                        'proteingym-
                        dms/PKN1_HUMAN_Tsuboyama_2023_1URF_indels',
                        'proteingym-dms/POLG_CXB3N_Mattenberger_2021',
                        'proteingym-dms/POLG_DEN26_Suphatrakul_2023',
                        'proteingym-dms/POLG_HCVJF_Qi_2014', 'proteingym-
                        dms/POLG_PESV_Tsuboyama_2023_2MXD', 'proteingym-
                        dms/POLG_PESV_Tsuboyama_2023_2MXD_indels',
                        'proteingym-dms/PPARG_HUMAN_Majithia_2016',
                        'proteingym-dms/PPM1D_HUMAN_Miller_2022', 'proteingym-
                        dms/PR40A_HUMAN_Tsuboyama_2023_1UZC', 'proteingym-
                        dms/PR40A_HUMAN_Tsuboyama_2023_1UZC_indels',
                        'proteingym-dms/PRKN_HUMAN_Clausen_2023', 'proteingym-
                        dms/PSAE_PICP2_Tsuboyama_2023_1PSE', 'proteingym-
                        dms/PSAE_PICP2_Tsuboyama_2023_1PSE_indels',
                        'proteingym-dms/PTEN_HUMAN_Matreyek_2021',
                        'proteingym-dms/PTEN_HUMAN_Mighell_2018', 'proteingym-
                        dms/PTEN_HUMAN_Mighell_2018_indels', 'proteingym-
                        dms/Q2N0S5_9HIV1_Haddox_2018', 'proteingym-
                        dms/Q53Z42_HUMAN_McShan_2019_binding-TAPBPR',
                        'proteingym-dms/Q53Z42_HUMAN_McShan_2019_expression',
                        'proteingym-dms/Q59976_STRSQ_Romero_2015',
                        'proteingym-dms/Q6WV12_9MAXI_Somermeyer_2022',
                        'proteingym-dms/Q837P4_ENTFA_Meier_2023', 'proteingym-
                        dms/Q837P5_ENTFA_Meier_2023', 'proteingym-
                        dms/Q8EG35_SHEON_Campbell_2022_indels', 'proteingym-
                        dms/Q8WTC7_9CNID_Somermeyer_2022', 'proteingym-
                        dms/R1AB_SARS2_Flynn_2022', 'proteingym-
                        dms/RAD_ANTMA_Tsuboyama_2023_2CJJ', 'proteingym-
                        dms/RAD_ANTMA_Tsuboyama_2023_2CJJ_indels',
                        'proteingym-dms/RAF1_HUMAN_Zinkus-Boltz_2019',
                        'proteingym-dms/RASH_HUMAN_Bandaru_2017', 'proteingym-
                        dms/RASK_HUMAN_Weng_2022_abundance', 'proteingym-
                        dms/RASK_HUMAN_Weng_2022_binding-DARPin_K55',
                        'proteingym-dms/RBP1_HUMAN_Tsuboyama_2023_2KWH',
                        'proteingym-dms/RCD1_ARATH_Tsuboyama_2023_5OAO',
                        'proteingym-
                        dms/RCD1_ARATH_Tsuboyama_2023_5OAO_indels',
                        'proteingym-dms/RCRO_LAMBD_Tsuboyama_2023_1ORC',
                        'proteingym-dms/RD23A_HUMAN_Tsuboyama_2023_1IFY',
                        'proteingym-
                        dms/RD23A_HUMAN_Tsuboyama_2023_1IFY_indels',
                        'proteingym-dms/RDRP_I33A0_Li_2023', 'proteingym-
                        dms/REV_HV1H2_Fernandes_2016', 'proteingym-
                        dms/RFAH_ECOLI_Tsuboyama_2023_2LCL', 'proteingym-
                        dms/RL20_AQUAE_Tsuboyama_2023_1GYZ', 'proteingym-
                        dms/RL40A_YEAST_Mavor_2016', 'proteingym-
                        dms/RL40A_YEAST_Roscoe_2013', 'proteingym-
                        dms/RL40A_YEAST_Roscoe_2014', 'proteingym-
                        dms/RNC_ECOLI_Weeks_2023', 'proteingym-
                        dms/RPC1_BP434_Tsuboyama_2023_1R69', 'proteingym-
                        dms/RPC1_BP434_Tsuboyama_2023_1R69_indels',
                        'proteingym-dms/RPC1_LAMBD_Li_2019_high-expression',
                        'proteingym-dms/RPC1_LAMBD_Li_2019_low-expression',
                        'proteingym-dms/RS15_GEOSE_Tsuboyama_2023_1A32',
                        'proteingym-
                        dms/RS15_GEOSE_Tsuboyama_2023_1A32_indels',
                        'proteingym-dms/S22A1_HUMAN_Yee_2023_abundance',
                        'proteingym-
                        dms/S22A1_HUMAN_Yee_2023_abundance_indels',
                        'proteingym-dms/S22A1_HUMAN_Yee_2023_activity',
                        'proteingym-dms/S22A1_HUMAN_Yee_2023_activity_indels',
                        'proteingym-dms/SAV1_MOUSE_Tsuboyama_2023_2YSB',
                        'proteingym-
                        dms/SAV1_MOUSE_Tsuboyama_2023_2YSB_indels',
                        'proteingym-dms/SBI_STAAM_Tsuboyama_2023_2JVG',
                        'proteingym-dms/SC6A4_HUMAN_Young_2021', 'proteingym-
                        dms/SCIN_STAAR_Tsuboyama_2023_2QFF', 'proteingym-
                        dms/SCN5A_HUMAN_Glazer_2019', 'proteingym-
                        dms/SDA_BACSU_Tsuboyama_2023_1PV0', 'proteingym-
                        dms/SDA_BACSU_Tsuboyama_2023_1PV0_indels',
                        'proteingym-dms/SERC_HUMAN_Xie_2023', 'proteingym-
                        dms/SHOC2_HUMAN_Kwon_2022', 'proteingym-
                        dms/SOX30_HUMAN_Tsuboyama_2023_7JJK', 'proteingym-
                        dms/SOX30_HUMAN_Tsuboyama_2023_7JJK_indels',
                        'proteingym-dms/SPA_STAAU_Tsuboyama_2023_1LP1',
                        'proteingym-dms/SPG1_STRSG_Olson_2014', 'proteingym-
                        dms/SPG1_STRSG_Wu_2016', 'proteingym-
                        dms/SPG2_STRSG_Tsuboyama_2023_5UBS', 'proteingym-
                        dms/SPG2_STRSG_Tsuboyama_2023_5UBS_indels',
                        'proteingym-dms/SPIKE_SARS2_Starr_2020_binding',
                        'proteingym-dms/SPIKE_SARS2_Starr_2020_expression',
                        'proteingym-dms/SPTN1_CHICK_Tsuboyama_2023_1TUD',
                        'proteingym-
                        dms/SPTN1_CHICK_Tsuboyama_2023_1TUD_indels',
                        'proteingym-dms/SQSTM_MOUSE_Tsuboyama_2023_2RRU',
                        'proteingym-
                        dms/SQSTM_MOUSE_Tsuboyama_2023_2RRU_indels',
                        'proteingym-dms/SR43C_ARATH_Tsuboyama_2023_2N88',
                        'proteingym-
                        dms/SR43C_ARATH_Tsuboyama_2023_2N88_indels',
                        'proteingym-dms/SRBS1_HUMAN_Tsuboyama_2023_2O2W',
                        'proteingym-
                        dms/SRBS1_HUMAN_Tsuboyama_2023_2O2W_indels',
                        'proteingym-dms/SRC_HUMAN_Ahler_2019', 'proteingym-
                        dms/SRC_HUMAN_Chakraborty_2023_binding-DAS_25uM',
                        'proteingym-dms/SRC_HUMAN_Nguyen_2022', 'proteingym-
                        dms/SUMO1_HUMAN_Weile_2017', 'proteingym-
                        dms/SYUA_HUMAN_Newberry_2020', 'proteingym-
                        dms/TADBP_HUMAN_Bolognesi_2019', 'proteingym-
                        dms/TAT_HV1BR_Fernandes_2016', 'proteingym-
                        dms/TCRG1_MOUSE_Tsuboyama_2023_1E0L', 'proteingym-
                        dms/TCRG1_MOUSE_Tsuboyama_2023_1E0L_indels',
                        'proteingym-dms/THO1_YEAST_Tsuboyama_2023_2WQG',
                        'proteingym-
                        dms/THO1_YEAST_Tsuboyama_2023_2WQG_indels',
                        'proteingym-dms/TNKS2_HUMAN_Tsuboyama_2023_5JRT',
                        'proteingym-
                        dms/TNKS2_HUMAN_Tsuboyama_2023_5JRT_indels',
                        'proteingym-dms/TPK1_HUMAN_Weile_2017', 'proteingym-
                        dms/TPMT_HUMAN_Matreyek_2018', 'proteingym-
                        dms/TPOR_HUMAN_Bridgford_2020', 'proteingym-
                        dms/TRPC_SACS2_Chan_2017', 'proteingym-
                        dms/TRPC_THEMA_Chan_2017', 'proteingym-
                        dms/UBC9_HUMAN_Weile_2017', 'proteingym-
                        dms/UBE4B_HUMAN_Tsuboyama_2023_3L1X', 'proteingym-
                        dms/UBE4B_HUMAN_Tsuboyama_2023_3L1X_indels',
                        'proteingym-dms/UBE4B_MOUSE_Starita_2013',
                        'proteingym-dms/UBR5_HUMAN_Tsuboyama_2023_1I2T',
                        'proteingym-
                        dms/UBR5_HUMAN_Tsuboyama_2023_1I2T_indels',
                        'proteingym-dms/VG08_BPP22_Tsuboyama_2023_2GP8',
                        'proteingym-
                        dms/VG08_BPP22_Tsuboyama_2023_2GP8_indels',
                        'proteingym-dms/VILI_CHICK_Tsuboyama_2023_1YU5',
                        'proteingym-
                        dms/VILI_CHICK_Tsuboyama_2023_1YU5_indels',
                        'proteingym-dms/VKOR1_HUMAN_Chiasson_2020_abundance',
                        'proteingym-dms/VKOR1_HUMAN_Chiasson_2020_activity',
                        'proteingym-dms/VRPI_BPT7_Tsuboyama_2023_2WNM',
                        'proteingym-dms/VRPI_BPT7_Tsuboyama_2023_2WNM_indels',
                        'proteingym-dms/YAIA_ECOLI_Tsuboyama_2023_2KVT',
                        'proteingym-dms/YAP1_HUMAN_Araya_2012', 'proteingym-
                        dms/YNZC_BACSU_Tsuboyama_2023_2JVD', 'proteingym-
                        dms/YNZC_BACSU_Tsuboyama_2023_2JVD_indels']
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
  --s3-bucket S3_BUCKET
                        S3 bucket name for uploading run artifacts (can also
                        set BIOMLBENCH_S3_BUCKET env var)
  --s3-prefix S3_PREFIX
                        S3 key prefix for uploaded artifacts (default: empty)
  --s3-no-compress      Disable compression when uploading to S3
  --disable-s3          Disable S3 uploads even if bucket is configured
  --incremental-s3      Enable incremental S3 uploads during execution (useful
                        for long benchmarks)
  --cpu-only            Use CPU-only container configuration (no GPU
                        requirements)
  --fast                Use fast container configuration (optimized entrypoint
                        for faster startup)

```

### `biomlbench grade`

```
usage: cli.py grade [-h] --submission SUBMISSION --output-dir OUTPUT_DIR
                    [--data-dir DATA_DIR] [--s3-bucket S3_BUCKET]
                    [--s3-prefix S3_PREFIX] [--s3-no-compress] [--disable-s3]

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
  --s3-bucket S3_BUCKET
                        S3 bucket name for uploading grading results (can also
                        set BIOMLBENCH_S3_BUCKET env var)
  --s3-prefix S3_PREFIX
                        S3 key prefix for uploaded artifacts (default: empty)
  --s3-no-compress      Disable compression when uploading to S3
  --disable-s3          Disable S3 uploads even if bucket is configured

```

### `biomlbench grade-sample`

```
usage: cli.py grade-sample [-h] [--data-dir DATA_DIR] submission task_id

positional arguments:
  submission           Path to the submission CSV file.
  task_id              ID of the task to grade in 'folder/task' format.
                       Examples: manual/caco2-wang

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
                        ['kaggle/histopathologic-cancer-detection',
                        'kaggle/hms-harmful-brain-activity-classification',
                        'kaggle/osic-pulmonary-fibrosis-progression',
                        'kaggle/ranzcr-clip-catheter-line-classification',
                        'kaggle/rsna-miccai-brain-tumor-radiogenomic-
                        classification', 'kaggle/stanford-covid-vaccine',
                        'kaggle/uw-madison-gi-tract-image-segmentation',
                        'kaggle/ventilator-pressure-prediction',
                        'kaggle_big/bms-molecular-translation',
                        'kaggle_big/rsna-2022-cervical-spine-fracture-
                        detection', 'kaggle_big/rsna-breast-cancer-detection',
                        'kaggle_big/siim-covid19-detection', 'kaggle_big/siim-
                        isic-melanoma-classification', 'kaggle_big/vinbigdata-
                        chest-xray-abnormalities-detection', 'manual/open-
                        problems-cell-cell-communication-ligand-target',
                        'manual/open-problems-label-projection', 'manual/open-
                        problems-predict-modality', 'manual/open-problems-
                        single-cell-perturbations', 'manual/open-problems-
                        spatially-variable-genes', 'polarishub/polaris-adme-
                        fang-hclint-1', 'polarishub/polaris-adme-fang-hppb-1',
                        'polarishub/polaris-adme-fang-perm-1',
                        'polarishub/polaris-adme-fang-r-1',
                        'polarishub/polaris-adme-fang-rclint-1',
                        'polarishub/polaris-adme-fang-rppb-1',
                        'polarishub/polaris-adme-fang-solu-1',
                        'polarishub/polaris-pkis1-egfr-wt-mut-c-1',
                        'polarishub/polaris-pkis1-egfr-wt-mut-r-1',
                        'polarishub/polaris-pkis1-kit-wt-mut-c-1',
                        'polarishub/polaris-pkis1-kit-wt-mut-r-1',
                        'polarishub/polaris-pkis1-ret-wt-mut-c-1',
                        'polarishub/polaris-pkis1-ret-wt-mut-r-1',
                        'polarishub/polaris-pkis2-egfr-wt-c-1',
                        'polarishub/polaris-pkis2-egfr-wt-r-1',
                        'polarishub/polaris-pkis2-kit-wt-c-1',
                        'polarishub/polaris-pkis2-kit-wt-r-1',
                        'polarishub/polaris-pkis2-lok-slk-c-1',
                        'polarishub/polaris-pkis2-lok-slk-r-1',
                        'polarishub/polaris-pkis2-ret-wt-c-1',
                        'polarishub/polaris-pkis2-ret-wt-r-1',
                        'polarishub/tdcommons-bbb-martins',
                        'polarishub/tdcommons-bioavailability-ma',
                        'polarishub/tdcommons-caco2-wang',
                        'polarishub/tdcommons-clearance-hepatocyte-az',
                        'polarishub/tdcommons-cyp2d6-substrate-carbonmangels',
                        'polarishub/tdcommons-herg', 'polarishub/tdcommons-
                        hia-hou', 'polarishub/tdcommons-lipophilicity-
                        astrazeneca', 'polarishub/tdcommons-pgp-broccatelli',
                        'polarishub-extended/adaptyv-bio-egfr-binders-binary-
                        cls-v0', 'polarishub-extended/biogen-adme-fang-hclint-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-hppb-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-perm-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-rclint-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-rppb-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-solu-
                        reg-v1', 'polarishub-extended/graphium-l1000-mcf7-v1',
                        'polarishub-extended/graphium-l1000-vcap-v1',
                        'polarishub-extended/graphium-pcba-1328-1564k-v1',
                        'polarishub-extended/graphium-qm9-v1', 'polarishub-
                        extended/graphium-tox21-v1', 'polarishub-
                        extended/graphium-zinc12k-v1', 'polarishub-
                        extended/novartis-adme-novartis-cyp3a4-cls',
                        'polarishub-extended/novartis-adme-novartis-
                        cyp3a4-reg', 'polarishub-extended/polaris-adme-fang-
                        hclint-1', 'polarishub-extended/polaris-adme-fang-
                        hppb-1', 'polarishub-extended/polaris-adme-fang-
                        perm-1', 'polarishub-extended/polaris-adme-fang-r-1',
                        'polarishub-extended/polaris-adme-fang-rclint-1',
                        'polarishub-extended/polaris-adme-fang-rppb-1',
                        'polarishub-extended/polaris-adme-fang-solu-1',
                        'polarishub-extended/polaris-hello-world-benchmark',
                        'polarishub-extended/polaris-
                        molprop-250k-leadlike-r-1', 'polarishub-
                        extended/polaris-molprop-250k-r-1', 'polarishub-
                        extended/polaris-molprop-250k-reg-v2', 'polarishub-
                        extended/polaris-molprop-leadlike-250k-reg-v2',
                        'polarishub-extended/polaris-pkis1-egfr-wt-mut-c-1',
                        'polarishub-extended/polaris-pkis1-egfr-wt-mut-r-1',
                        'polarishub-extended/polaris-pkis1-kit-wt-mut-c-1',
                        'polarishub-extended/polaris-pkis1-kit-wt-mut-r-1',
                        'polarishub-extended/polaris-pkis1-ret-wt-mut-c-1',
                        'polarishub-extended/polaris-pkis1-ret-wt-mut-r-1',
                        'polarishub-extended/polaris-pkis2-egfr-wt-c-1',
                        'polarishub-extended/polaris-pkis2-egfr-wt-r-1',
                        'polarishub-extended/polaris-pkis2-egfr-wt-reg-v2',
                        'polarishub-extended/polaris-pkis2-kit-wt-c-1',
                        'polarishub-extended/polaris-pkis2-kit-wt-cls-v2',
                        'polarishub-extended/polaris-pkis2-kit-wt-r-1',
                        'polarishub-extended/polaris-pkis2-kit-wt-reg-v2',
                        'polarishub-extended/polaris-pkis2-lok-slk-c-1',
                        'polarishub-extended/polaris-pkis2-lok-slk-cls-v2',
                        'polarishub-extended/polaris-pkis2-lok-slk-r-1',
                        'polarishub-extended/polaris-pkis2-ret-wt-c-1',
                        'polarishub-extended/polaris-pkis2-ret-wt-cls-v2',
                        'polarishub-extended/polaris-pkis2-ret-wt-r-1',
                        'polarishub-extended/polaris-pkis2-ret-wt-reg-v2',
                        'polarishub-extended/polaris-posebusters-v1',
                        'polarishub-extended/recursion-rxrx-compound-gene-
                        activity-benchmark', 'polarishub-extended/tdcommons-
                        ames', 'polarishub-extended/tdcommons-bbb-martins',
                        'polarishub-extended/tdcommons-bioavailability-ma',
                        'polarishub-extended/tdcommons-caco2-wang',
                        'polarishub-extended/tdcommons-clearance-hepatocyte-
                        az', 'polarishub-extended/tdcommons-clearance-
                        microsome-az', 'polarishub-extended/tdcommons-
                        cyp2c9-substrate-carbonmangels', 'polarishub-
                        extended/tdcommons-cyp2c9-veith', 'polarishub-
                        extended/tdcommons-cyp2d6-substrate-carbonmangels',
                        'polarishub-extended/tdcommons-cyp2d6-veith',
                        'polarishub-extended/tdcommons-cyp3a4-substrate-
                        carbonmangels', 'polarishub-extended/tdcommons-
                        cyp3a4-veith', 'polarishub-extended/tdcommons-dili',
                        'polarishub-extended/tdcommons-half-life-obach',
                        'polarishub-extended/tdcommons-herg', 'polarishub-
                        extended/tdcommons-hia-hou', 'polarishub-
                        extended/tdcommons-ld50-zhu', 'polarishub-
                        extended/tdcommons-lipophilicity-astrazeneca',
                        'polarishub-extended/tdcommons-pgp-broccatelli',
                        'polarishub-extended/tdcommons-ppbr-az', 'polarishub-
                        extended/tdcommons-solubility-aqsoldb', 'polarishub-
                        extended/tdcommons-vdss-lombardo', 'polarishub-
                        extended/vishrut64-cho-dna-expression-prediction-
                        dataset-task', 'polarishub-extended/vishrut64-rna-
                        expression-prediction-dataset-task', 'proteingym-
                        dms/A0A140D2T1_ZIKV_Sourisseau_2019', 'proteingym-
                        dms/A0A192B1T2_9HIV1_Haddox_2018', 'proteingym-
                        dms/A0A1I9GEU1_NEIME_Kennouche_2019', 'proteingym-
                        dms/A0A247D711_LISMN_Stadelmann_2021', 'proteingym-
                        dms/A0A2Z5U3Z0_9INFA_Doud_2016', 'proteingym-
                        dms/A0A2Z5U3Z0_9INFA_Wu_2014', 'proteingym-
                        dms/A4D664_9INFA_Soh_2019', 'proteingym-
                        dms/A4GRB6_PSEAI_Chen_2020', 'proteingym-
                        dms/A4_HUMAN_Seuma_2022', 'proteingym-
                        dms/A4_HUMAN_Seuma_2022_indels', 'proteingym-
                        dms/AACC1_PSEAI_Dandage_2018', 'proteingym-
                        dms/ACE2_HUMAN_Chan_2020', 'proteingym-
                        dms/ADRB2_HUMAN_Jones_2020', 'proteingym-
                        dms/AICDA_HUMAN_Gajula_2014_3cycles', 'proteingym-
                        dms/AMFR_HUMAN_Tsuboyama_2023_4G3O', 'proteingym-
                        dms/AMFR_HUMAN_Tsuboyama_2023_4G3O_indels',
                        'proteingym-dms/AMIE_PSEAE_Wrenbeck_2017',
                        'proteingym-dms/ANCSZ_Hobbs_2022', 'proteingym-
                        dms/ARGR_ECOLI_Tsuboyama_2023_1AOY', 'proteingym-
                        dms/ARGR_ECOLI_Tsuboyama_2023_1AOY_indels',
                        'proteingym-dms/B1LPA6_ECOSM_Russ_2020_indels',
                        'proteingym-dms/B2L11_HUMAN_Dutta_2010_binding-Mcl-1',
                        'proteingym-dms/BBC1_YEAST_Tsuboyama_2023_1TG0',
                        'proteingym-
                        dms/BBC1_YEAST_Tsuboyama_2023_1TG0_indels',
                        'proteingym-dms/BCHB_CHLTE_Tsuboyama_2023_2KRU',
                        'proteingym-
                        dms/BCHB_CHLTE_Tsuboyama_2023_2KRU_indels',
                        'proteingym-dms/BLAT_ECOLX_Deng_2012', 'proteingym-
                        dms/BLAT_ECOLX_Firnberg_2014', 'proteingym-
                        dms/BLAT_ECOLX_Gonzalez_2019_indels', 'proteingym-
                        dms/BLAT_ECOLX_Jacquier_2013', 'proteingym-
                        dms/BLAT_ECOLX_Stiffler_2015', 'proteingym-
                        dms/BRCA1_HUMAN_Findlay_2018', 'proteingym-
                        dms/BRCA2_HUMAN_Erwood_2022_HEK293T', 'proteingym-
                        dms/C6KNH7_9INFA_Lee_2018', 'proteingym-
                        dms/CALM1_HUMAN_Weile_2017', 'proteingym-
                        dms/CAPSD_AAV2S_Sinai_2021', 'proteingym-
                        dms/CAPSD_AAV2S_Sinai_2021_designed_indels',
                        'proteingym-
                        dms/CAPSD_AAV2S_Sinai_2021_library_indels',
                        'proteingym-dms/CAR11_HUMAN_Meitlis_2020_gof',
                        'proteingym-dms/CAR11_HUMAN_Meitlis_2020_lof',
                        'proteingym-dms/CAS9_STRP1_Spencer_2017_positive',
                        'proteingym-dms/CASP3_HUMAN_Roychowdhury_2020',
                        'proteingym-dms/CASP7_HUMAN_Roychowdhury_2020',
                        'proteingym-dms/CATR_CHLRE_Tsuboyama_2023_2AMI',
                        'proteingym-
                        dms/CATR_CHLRE_Tsuboyama_2023_2AMI_indels',
                        'proteingym-dms/CBPA2_HUMAN_Tsuboyama_2023_1O6X',
                        'proteingym-
                        dms/CBPA2_HUMAN_Tsuboyama_2023_1O6X_indels',
                        'proteingym-dms/CBS_HUMAN_Sun_2020', 'proteingym-
                        dms/CBX4_HUMAN_Tsuboyama_2023_2K28', 'proteingym-
                        dms/CBX4_HUMAN_Tsuboyama_2023_2K28_indels',
                        'proteingym-dms/CCDB_ECOLI_Adkar_2012', 'proteingym-
                        dms/CCDB_ECOLI_Tripathi_2016', 'proteingym-
                        dms/CCR5_HUMAN_Gill_2023', 'proteingym-
                        dms/CD19_HUMAN_Klesmith_2019_FMC_singles',
                        'proteingym-dms/CP2C9_HUMAN_Amorosi_2021_abundance',
                        'proteingym-dms/CP2C9_HUMAN_Amorosi_2021_activity',
                        'proteingym-dms/CSN4_MOUSE_Tsuboyama_2023_1UFM',
                        'proteingym-
                        dms/CSN4_MOUSE_Tsuboyama_2023_1UFM_indels',
                        'proteingym-dms/CUE1_YEAST_Tsuboyama_2023_2MYX',
                        'proteingym-
                        dms/CUE1_YEAST_Tsuboyama_2023_2MYX_indels',
                        'proteingym-dms/D7PM05_CLYGR_Somermeyer_2022',
                        'proteingym-dms/DLG4_HUMAN_Faure_2021', 'proteingym-
                        dms/DLG4_RAT_McLaughlin_2012', 'proteingym-
                        dms/DN7A_SACS2_Tsuboyama_2023_1JIC', 'proteingym-
                        dms/DN7A_SACS2_Tsuboyama_2023_1JIC_indels',
                        'proteingym-dms/DNJA1_HUMAN_Tsuboyama_2023_2LO1',
                        'proteingym-
                        dms/DNJA1_HUMAN_Tsuboyama_2023_2LO1_indels',
                        'proteingym-dms/DOCK1_MOUSE_Tsuboyama_2023_2M0Y',
                        'proteingym-
                        dms/DOCK1_MOUSE_Tsuboyama_2023_2M0Y_indels',
                        'proteingym-dms/DYR_ECOLI_Nguyen_2023', 'proteingym-
                        dms/DYR_ECOLI_Thompson_2019', 'proteingym-
                        dms/ENVZ_ECOLI_Ghose_2023', 'proteingym-
                        dms/ENV_HV1B9_DuenasDecamp_2016', 'proteingym-
                        dms/ENV_HV1BR_Haddox_2016', 'proteingym-
                        dms/EPHB2_HUMAN_Tsuboyama_2023_1F0M', 'proteingym-
                        dms/EPHB2_HUMAN_Tsuboyama_2023_1F0M_indels',
                        'proteingym-dms/ERBB2_HUMAN_Elazar_2016', 'proteingym-
                        dms/ESTA_BACSU_Nutschel_2020', 'proteingym-
                        dms/F7YBW8_MESOW_Aakre_2015', 'proteingym-
                        dms/F7YBW8_MESOW_Ding_2023', 'proteingym-
                        dms/FECA_ECOLI_Tsuboyama_2023_2D1U', 'proteingym-
                        dms/FECA_ECOLI_Tsuboyama_2023_2D1U_indels',
                        'proteingym-dms/FKBP3_HUMAN_Tsuboyama_2023_2KFV',
                        'proteingym-dms/GAL4_YEAST_Kitzman_2015', 'proteingym-
                        dms/GCN4_YEAST_Staller_2018', 'proteingym-
                        dms/GDIA_HUMAN_Silverstein_2021', 'proteingym-
                        dms/GFP_AEQVI_Sarkisyan_2016', 'proteingym-
                        dms/GLPA_HUMAN_Elazar_2016', 'proteingym-
                        dms/GRB2_HUMAN_Faure_2021', 'proteingym-
                        dms/HCP_LAMBD_Tsuboyama_2023_2L6Q', 'proteingym-
                        dms/HCP_LAMBD_Tsuboyama_2023_2L6Q_indels',
                        'proteingym-dms/HECD1_HUMAN_Tsuboyama_2023_3DKM',
                        'proteingym-
                        dms/HECD1_HUMAN_Tsuboyama_2023_3DKM_indels',
                        'proteingym-dms/HEM3_HUMAN_Loggerenberg_2023',
                        'proteingym-dms/HIS7_YEAST_Pokusaeva_2019',
                        'proteingym-dms/HIS7_YEAST_Pokusaeva_2019_indels',
                        'proteingym-dms/HMDH_HUMAN_Jiang_2019', 'proteingym-
                        dms/HSP82_YEAST_Cote-Hammarlof_2020_growth-H2O2',
                        'proteingym-dms/HSP82_YEAST_Flynn_2019', 'proteingym-
                        dms/HSP82_YEAST_Mishra_2016', 'proteingym-
                        dms/HXK4_HUMAN_Gersing_2022_activity', 'proteingym-
                        dms/HXK4_HUMAN_Gersing_2023_abundance', 'proteingym-
                        dms/I6TAH8_I68A0_Doud_2015', 'proteingym-
                        dms/IF1_ECOLI_Kelsic_2016', 'proteingym-
                        dms/ILF3_HUMAN_Tsuboyama_2023_2L33', 'proteingym-
                        dms/ILF3_HUMAN_Tsuboyama_2023_2L33_indels',
                        'proteingym-dms/ISDH_STAAW_Tsuboyama_2023_2LHR',
                        'proteingym-dms/KCNE1_HUMAN_Muhammad_2023_expression',
                        'proteingym-dms/KCNE1_HUMAN_Muhammad_2023_function',
                        'proteingym-dms/KCNH2_HUMAN_Kozek_2020', 'proteingym-
                        dms/KCNJ2_MOUSE_Coyote-Maestas_2022_function',
                        'proteingym-dms/KCNJ2_MOUSE_Coyote-
                        Maestas_2022_surface', 'proteingym-
                        dms/KCNJ2_MOUSE_Macdonald_2022_indels', 'proteingym-
                        dms/KKA2_KLEPN_Melnikov_2014', 'proteingym-
                        dms/LGK_LIPST_Klesmith_2015', 'proteingym-
                        dms/LYAM1_HUMAN_Elazar_2016', 'proteingym-
                        dms/MAFG_MOUSE_Tsuboyama_2023_1K1V', 'proteingym-
                        dms/MAFG_MOUSE_Tsuboyama_2023_1K1V_indels',
                        'proteingym-dms/MBD11_ARATH_Tsuboyama_2023_6ACV',
                        'proteingym-
                        dms/MBD11_ARATH_Tsuboyama_2023_6ACV_indels',
                        'proteingym-dms/MET_HUMAN_Estevam_2023', 'proteingym-
                        dms/MK01_HUMAN_Brenan_2016', 'proteingym-
                        dms/MLAC_ECOLI_MacRae_2023', 'proteingym-
                        dms/MSH2_HUMAN_Jia_2020', 'proteingym-
                        dms/MTH3_HAEAE_RockahShmuel_2015', 'proteingym-
                        dms/MTHR_HUMAN_Weile_2021', 'proteingym-
                        dms/MYO3_YEAST_Tsuboyama_2023_2BTT', 'proteingym-
                        dms/MYO3_YEAST_Tsuboyama_2023_2BTT_indels',
                        'proteingym-dms/NCAP_I34A1_Doud_2015', 'proteingym-
                        dms/NKX31_HUMAN_Tsuboyama_2023_2L9R', 'proteingym-
                        dms/NKX31_HUMAN_Tsuboyama_2023_2L9R_indels',
                        'proteingym-dms/NPC1_HUMAN_Erwood_2022_HEK293T',
                        'proteingym-dms/NPC1_HUMAN_Erwood_2022_RPE1',
                        'proteingym-dms/NRAM_I33A0_Jiang_2016', 'proteingym-
                        dms/NUD15_HUMAN_Suiter_2020', 'proteingym-
                        dms/NUSA_ECOLI_Tsuboyama_2023_1WCL', 'proteingym-
                        dms/NUSA_ECOLI_Tsuboyama_2023_1WCL_indels',
                        'proteingym-dms/NUSG_MYCTU_Tsuboyama_2023_2MI6',
                        'proteingym-
                        dms/NUSG_MYCTU_Tsuboyama_2023_2MI6_indels',
                        'proteingym-dms/OBSCN_HUMAN_Tsuboyama_2023_1V1C',
                        'proteingym-
                        dms/OBSCN_HUMAN_Tsuboyama_2023_1V1C_indels',
                        'proteingym-dms/ODP2_GEOSE_Tsuboyama_2023_1W4G',
                        'proteingym-
                        dms/ODP2_GEOSE_Tsuboyama_2023_1W4G_indels',
                        'proteingym-dms/OPSD_HUMAN_Wan_2019', 'proteingym-
                        dms/OTC_HUMAN_Lo_2023', 'proteingym-
                        dms/OTU7A_HUMAN_Tsuboyama_2023_2L2D', 'proteingym-
                        dms/OTU7A_HUMAN_Tsuboyama_2023_2L2D_indels',
                        'proteingym-dms/OXDA_RHOTO_Vanella_2023_activity',
                        'proteingym-dms/OXDA_RHOTO_Vanella_2023_expression',
                        'proteingym-
                        dms/P53_HUMAN_Giacomelli_2018_Null_Etoposide',
                        'proteingym-
                        dms/P53_HUMAN_Giacomelli_2018_Null_Nutlin',
                        'proteingym-dms/P53_HUMAN_Giacomelli_2018_WT_Nutlin',
                        'proteingym-dms/P53_HUMAN_Kotler_2018', 'proteingym-
                        dms/P53_HUMAN_Kotler_2018_indels', 'proteingym-
                        dms/P84126_THETH_Chan_2017', 'proteingym-
                        dms/PABP_YEAST_Melamed_2013', 'proteingym-
                        dms/PAI1_HUMAN_Huttinger_2021', 'proteingym-
                        dms/PA_I34A1_Wu_2015', 'proteingym-
                        dms/PHOT_CHLRE_Chen_2023', 'proteingym-
                        dms/PIN1_HUMAN_Tsuboyama_2023_1I6C', 'proteingym-
                        dms/PIN1_HUMAN_Tsuboyama_2023_1I6C_indels',
                        'proteingym-dms/PITX2_HUMAN_Tsuboyama_2023_2L7M',
                        'proteingym-
                        dms/PITX2_HUMAN_Tsuboyama_2023_2L7M_indels',
                        'proteingym-dms/PKN1_HUMAN_Tsuboyama_2023_1URF',
                        'proteingym-
                        dms/PKN1_HUMAN_Tsuboyama_2023_1URF_indels',
                        'proteingym-dms/POLG_CXB3N_Mattenberger_2021',
                        'proteingym-dms/POLG_DEN26_Suphatrakul_2023',
                        'proteingym-dms/POLG_HCVJF_Qi_2014', 'proteingym-
                        dms/POLG_PESV_Tsuboyama_2023_2MXD', 'proteingym-
                        dms/POLG_PESV_Tsuboyama_2023_2MXD_indels',
                        'proteingym-dms/PPARG_HUMAN_Majithia_2016',
                        'proteingym-dms/PPM1D_HUMAN_Miller_2022', 'proteingym-
                        dms/PR40A_HUMAN_Tsuboyama_2023_1UZC', 'proteingym-
                        dms/PR40A_HUMAN_Tsuboyama_2023_1UZC_indels',
                        'proteingym-dms/PRKN_HUMAN_Clausen_2023', 'proteingym-
                        dms/PSAE_PICP2_Tsuboyama_2023_1PSE', 'proteingym-
                        dms/PSAE_PICP2_Tsuboyama_2023_1PSE_indels',
                        'proteingym-dms/PTEN_HUMAN_Matreyek_2021',
                        'proteingym-dms/PTEN_HUMAN_Mighell_2018', 'proteingym-
                        dms/PTEN_HUMAN_Mighell_2018_indels', 'proteingym-
                        dms/Q2N0S5_9HIV1_Haddox_2018', 'proteingym-
                        dms/Q53Z42_HUMAN_McShan_2019_binding-TAPBPR',
                        'proteingym-dms/Q53Z42_HUMAN_McShan_2019_expression',
                        'proteingym-dms/Q59976_STRSQ_Romero_2015',
                        'proteingym-dms/Q6WV12_9MAXI_Somermeyer_2022',
                        'proteingym-dms/Q837P4_ENTFA_Meier_2023', 'proteingym-
                        dms/Q837P5_ENTFA_Meier_2023', 'proteingym-
                        dms/Q8EG35_SHEON_Campbell_2022_indels', 'proteingym-
                        dms/Q8WTC7_9CNID_Somermeyer_2022', 'proteingym-
                        dms/R1AB_SARS2_Flynn_2022', 'proteingym-
                        dms/RAD_ANTMA_Tsuboyama_2023_2CJJ', 'proteingym-
                        dms/RAD_ANTMA_Tsuboyama_2023_2CJJ_indels',
                        'proteingym-dms/RAF1_HUMAN_Zinkus-Boltz_2019',
                        'proteingym-dms/RASH_HUMAN_Bandaru_2017', 'proteingym-
                        dms/RASK_HUMAN_Weng_2022_abundance', 'proteingym-
                        dms/RASK_HUMAN_Weng_2022_binding-DARPin_K55',
                        'proteingym-dms/RBP1_HUMAN_Tsuboyama_2023_2KWH',
                        'proteingym-dms/RCD1_ARATH_Tsuboyama_2023_5OAO',
                        'proteingym-
                        dms/RCD1_ARATH_Tsuboyama_2023_5OAO_indels',
                        'proteingym-dms/RCRO_LAMBD_Tsuboyama_2023_1ORC',
                        'proteingym-dms/RD23A_HUMAN_Tsuboyama_2023_1IFY',
                        'proteingym-
                        dms/RD23A_HUMAN_Tsuboyama_2023_1IFY_indels',
                        'proteingym-dms/RDRP_I33A0_Li_2023', 'proteingym-
                        dms/REV_HV1H2_Fernandes_2016', 'proteingym-
                        dms/RFAH_ECOLI_Tsuboyama_2023_2LCL', 'proteingym-
                        dms/RL20_AQUAE_Tsuboyama_2023_1GYZ', 'proteingym-
                        dms/RL40A_YEAST_Mavor_2016', 'proteingym-
                        dms/RL40A_YEAST_Roscoe_2013', 'proteingym-
                        dms/RL40A_YEAST_Roscoe_2014', 'proteingym-
                        dms/RNC_ECOLI_Weeks_2023', 'proteingym-
                        dms/RPC1_BP434_Tsuboyama_2023_1R69', 'proteingym-
                        dms/RPC1_BP434_Tsuboyama_2023_1R69_indels',
                        'proteingym-dms/RPC1_LAMBD_Li_2019_high-expression',
                        'proteingym-dms/RPC1_LAMBD_Li_2019_low-expression',
                        'proteingym-dms/RS15_GEOSE_Tsuboyama_2023_1A32',
                        'proteingym-
                        dms/RS15_GEOSE_Tsuboyama_2023_1A32_indels',
                        'proteingym-dms/S22A1_HUMAN_Yee_2023_abundance',
                        'proteingym-
                        dms/S22A1_HUMAN_Yee_2023_abundance_indels',
                        'proteingym-dms/S22A1_HUMAN_Yee_2023_activity',
                        'proteingym-dms/S22A1_HUMAN_Yee_2023_activity_indels',
                        'proteingym-dms/SAV1_MOUSE_Tsuboyama_2023_2YSB',
                        'proteingym-
                        dms/SAV1_MOUSE_Tsuboyama_2023_2YSB_indels',
                        'proteingym-dms/SBI_STAAM_Tsuboyama_2023_2JVG',
                        'proteingym-dms/SC6A4_HUMAN_Young_2021', 'proteingym-
                        dms/SCIN_STAAR_Tsuboyama_2023_2QFF', 'proteingym-
                        dms/SCN5A_HUMAN_Glazer_2019', 'proteingym-
                        dms/SDA_BACSU_Tsuboyama_2023_1PV0', 'proteingym-
                        dms/SDA_BACSU_Tsuboyama_2023_1PV0_indels',
                        'proteingym-dms/SERC_HUMAN_Xie_2023', 'proteingym-
                        dms/SHOC2_HUMAN_Kwon_2022', 'proteingym-
                        dms/SOX30_HUMAN_Tsuboyama_2023_7JJK', 'proteingym-
                        dms/SOX30_HUMAN_Tsuboyama_2023_7JJK_indels',
                        'proteingym-dms/SPA_STAAU_Tsuboyama_2023_1LP1',
                        'proteingym-dms/SPG1_STRSG_Olson_2014', 'proteingym-
                        dms/SPG1_STRSG_Wu_2016', 'proteingym-
                        dms/SPG2_STRSG_Tsuboyama_2023_5UBS', 'proteingym-
                        dms/SPG2_STRSG_Tsuboyama_2023_5UBS_indels',
                        'proteingym-dms/SPIKE_SARS2_Starr_2020_binding',
                        'proteingym-dms/SPIKE_SARS2_Starr_2020_expression',
                        'proteingym-dms/SPTN1_CHICK_Tsuboyama_2023_1TUD',
                        'proteingym-
                        dms/SPTN1_CHICK_Tsuboyama_2023_1TUD_indels',
                        'proteingym-dms/SQSTM_MOUSE_Tsuboyama_2023_2RRU',
                        'proteingym-
                        dms/SQSTM_MOUSE_Tsuboyama_2023_2RRU_indels',
                        'proteingym-dms/SR43C_ARATH_Tsuboyama_2023_2N88',
                        'proteingym-
                        dms/SR43C_ARATH_Tsuboyama_2023_2N88_indels',
                        'proteingym-dms/SRBS1_HUMAN_Tsuboyama_2023_2O2W',
                        'proteingym-
                        dms/SRBS1_HUMAN_Tsuboyama_2023_2O2W_indels',
                        'proteingym-dms/SRC_HUMAN_Ahler_2019', 'proteingym-
                        dms/SRC_HUMAN_Chakraborty_2023_binding-DAS_25uM',
                        'proteingym-dms/SRC_HUMAN_Nguyen_2022', 'proteingym-
                        dms/SUMO1_HUMAN_Weile_2017', 'proteingym-
                        dms/SYUA_HUMAN_Newberry_2020', 'proteingym-
                        dms/TADBP_HUMAN_Bolognesi_2019', 'proteingym-
                        dms/TAT_HV1BR_Fernandes_2016', 'proteingym-
                        dms/TCRG1_MOUSE_Tsuboyama_2023_1E0L', 'proteingym-
                        dms/TCRG1_MOUSE_Tsuboyama_2023_1E0L_indels',
                        'proteingym-dms/THO1_YEAST_Tsuboyama_2023_2WQG',
                        'proteingym-
                        dms/THO1_YEAST_Tsuboyama_2023_2WQG_indels',
                        'proteingym-dms/TNKS2_HUMAN_Tsuboyama_2023_5JRT',
                        'proteingym-
                        dms/TNKS2_HUMAN_Tsuboyama_2023_5JRT_indels',
                        'proteingym-dms/TPK1_HUMAN_Weile_2017', 'proteingym-
                        dms/TPMT_HUMAN_Matreyek_2018', 'proteingym-
                        dms/TPOR_HUMAN_Bridgford_2020', 'proteingym-
                        dms/TRPC_SACS2_Chan_2017', 'proteingym-
                        dms/TRPC_THEMA_Chan_2017', 'proteingym-
                        dms/UBC9_HUMAN_Weile_2017', 'proteingym-
                        dms/UBE4B_HUMAN_Tsuboyama_2023_3L1X', 'proteingym-
                        dms/UBE4B_HUMAN_Tsuboyama_2023_3L1X_indels',
                        'proteingym-dms/UBE4B_MOUSE_Starita_2013',
                        'proteingym-dms/UBR5_HUMAN_Tsuboyama_2023_1I2T',
                        'proteingym-
                        dms/UBR5_HUMAN_Tsuboyama_2023_1I2T_indels',
                        'proteingym-dms/VG08_BPP22_Tsuboyama_2023_2GP8',
                        'proteingym-
                        dms/VG08_BPP22_Tsuboyama_2023_2GP8_indels',
                        'proteingym-dms/VILI_CHICK_Tsuboyama_2023_1YU5',
                        'proteingym-
                        dms/VILI_CHICK_Tsuboyama_2023_1YU5_indels',
                        'proteingym-dms/VKOR1_HUMAN_Chiasson_2020_abundance',
                        'proteingym-dms/VKOR1_HUMAN_Chiasson_2020_activity',
                        'proteingym-dms/VRPI_BPT7_Tsuboyama_2023_2WNM',
                        'proteingym-dms/VRPI_BPT7_Tsuboyama_2023_2WNM_indels',
                        'proteingym-dms/YAIA_ECOLI_Tsuboyama_2023_2KVT',
                        'proteingym-dms/YAP1_HUMAN_Araya_2012', 'proteingym-
                        dms/YNZC_BACSU_Tsuboyama_2023_2JVD', 'proteingym-
                        dms/YNZC_BACSU_Tsuboyama_2023_2JVD_indels']

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
                        Valid options: ['kaggle/histopathologic-cancer-
                        detection', 'kaggle/hms-harmful-brain-activity-
                        classification', 'kaggle/osic-pulmonary-fibrosis-
                        progression', 'kaggle/ranzcr-clip-catheter-line-
                        classification', 'kaggle/rsna-miccai-brain-tumor-
                        radiogenomic-classification', 'kaggle/stanford-covid-
                        vaccine', 'kaggle/uw-madison-gi-tract-image-
                        segmentation', 'kaggle/ventilator-pressure-
                        prediction', 'kaggle_big/bms-molecular-translation',
                        'kaggle_big/rsna-2022-cervical-spine-fracture-
                        detection', 'kaggle_big/rsna-breast-cancer-detection',
                        'kaggle_big/siim-covid19-detection', 'kaggle_big/siim-
                        isic-melanoma-classification', 'kaggle_big/vinbigdata-
                        chest-xray-abnormalities-detection', 'manual/open-
                        problems-cell-cell-communication-ligand-target',
                        'manual/open-problems-label-projection', 'manual/open-
                        problems-predict-modality', 'manual/open-problems-
                        single-cell-perturbations', 'manual/open-problems-
                        spatially-variable-genes', 'polarishub/polaris-adme-
                        fang-hclint-1', 'polarishub/polaris-adme-fang-hppb-1',
                        'polarishub/polaris-adme-fang-perm-1',
                        'polarishub/polaris-adme-fang-r-1',
                        'polarishub/polaris-adme-fang-rclint-1',
                        'polarishub/polaris-adme-fang-rppb-1',
                        'polarishub/polaris-adme-fang-solu-1',
                        'polarishub/polaris-pkis1-egfr-wt-mut-c-1',
                        'polarishub/polaris-pkis1-egfr-wt-mut-r-1',
                        'polarishub/polaris-pkis1-kit-wt-mut-c-1',
                        'polarishub/polaris-pkis1-kit-wt-mut-r-1',
                        'polarishub/polaris-pkis1-ret-wt-mut-c-1',
                        'polarishub/polaris-pkis1-ret-wt-mut-r-1',
                        'polarishub/polaris-pkis2-egfr-wt-c-1',
                        'polarishub/polaris-pkis2-egfr-wt-r-1',
                        'polarishub/polaris-pkis2-kit-wt-c-1',
                        'polarishub/polaris-pkis2-kit-wt-r-1',
                        'polarishub/polaris-pkis2-lok-slk-c-1',
                        'polarishub/polaris-pkis2-lok-slk-r-1',
                        'polarishub/polaris-pkis2-ret-wt-c-1',
                        'polarishub/polaris-pkis2-ret-wt-r-1',
                        'polarishub/tdcommons-bbb-martins',
                        'polarishub/tdcommons-bioavailability-ma',
                        'polarishub/tdcommons-caco2-wang',
                        'polarishub/tdcommons-clearance-hepatocyte-az',
                        'polarishub/tdcommons-cyp2d6-substrate-carbonmangels',
                        'polarishub/tdcommons-herg', 'polarishub/tdcommons-
                        hia-hou', 'polarishub/tdcommons-lipophilicity-
                        astrazeneca', 'polarishub/tdcommons-pgp-broccatelli',
                        'polarishub-extended/adaptyv-bio-egfr-binders-binary-
                        cls-v0', 'polarishub-extended/biogen-adme-fang-hclint-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-hppb-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-perm-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-rclint-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-rppb-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-solu-
                        reg-v1', 'polarishub-extended/graphium-l1000-mcf7-v1',
                        'polarishub-extended/graphium-l1000-vcap-v1',
                        'polarishub-extended/graphium-pcba-1328-1564k-v1',
                        'polarishub-extended/graphium-qm9-v1', 'polarishub-
                        extended/graphium-tox21-v1', 'polarishub-
                        extended/graphium-zinc12k-v1', 'polarishub-
                        extended/novartis-adme-novartis-cyp3a4-cls',
                        'polarishub-extended/novartis-adme-novartis-
                        cyp3a4-reg', 'polarishub-extended/polaris-adme-fang-
                        hclint-1', 'polarishub-extended/polaris-adme-fang-
                        hppb-1', 'polarishub-extended/polaris-adme-fang-
                        perm-1', 'polarishub-extended/polaris-adme-fang-r-1',
                        'polarishub-extended/polaris-adme-fang-rclint-1',
                        'polarishub-extended/polaris-adme-fang-rppb-1',
                        'polarishub-extended/polaris-adme-fang-solu-1',
                        'polarishub-extended/polaris-hello-world-benchmark',
                        'polarishub-extended/polaris-
                        molprop-250k-leadlike-r-1', 'polarishub-
                        extended/polaris-molprop-250k-r-1', 'polarishub-
                        extended/polaris-molprop-250k-reg-v2', 'polarishub-
                        extended/polaris-molprop-leadlike-250k-reg-v2',
                        'polarishub-extended/polaris-pkis1-egfr-wt-mut-c-1',
                        'polarishub-extended/polaris-pkis1-egfr-wt-mut-r-1',
                        'polarishub-extended/polaris-pkis1-kit-wt-mut-c-1',
                        'polarishub-extended/polaris-pkis1-kit-wt-mut-r-1',
                        'polarishub-extended/polaris-pkis1-ret-wt-mut-c-1',
                        'polarishub-extended/polaris-pkis1-ret-wt-mut-r-1',
                        'polarishub-extended/polaris-pkis2-egfr-wt-c-1',
                        'polarishub-extended/polaris-pkis2-egfr-wt-r-1',
                        'polarishub-extended/polaris-pkis2-egfr-wt-reg-v2',
                        'polarishub-extended/polaris-pkis2-kit-wt-c-1',
                        'polarishub-extended/polaris-pkis2-kit-wt-cls-v2',
                        'polarishub-extended/polaris-pkis2-kit-wt-r-1',
                        'polarishub-extended/polaris-pkis2-kit-wt-reg-v2',
                        'polarishub-extended/polaris-pkis2-lok-slk-c-1',
                        'polarishub-extended/polaris-pkis2-lok-slk-cls-v2',
                        'polarishub-extended/polaris-pkis2-lok-slk-r-1',
                        'polarishub-extended/polaris-pkis2-ret-wt-c-1',
                        'polarishub-extended/polaris-pkis2-ret-wt-cls-v2',
                        'polarishub-extended/polaris-pkis2-ret-wt-r-1',
                        'polarishub-extended/polaris-pkis2-ret-wt-reg-v2',
                        'polarishub-extended/polaris-posebusters-v1',
                        'polarishub-extended/recursion-rxrx-compound-gene-
                        activity-benchmark', 'polarishub-extended/tdcommons-
                        ames', 'polarishub-extended/tdcommons-bbb-martins',
                        'polarishub-extended/tdcommons-bioavailability-ma',
                        'polarishub-extended/tdcommons-caco2-wang',
                        'polarishub-extended/tdcommons-clearance-hepatocyte-
                        az', 'polarishub-extended/tdcommons-clearance-
                        microsome-az', 'polarishub-extended/tdcommons-
                        cyp2c9-substrate-carbonmangels', 'polarishub-
                        extended/tdcommons-cyp2c9-veith', 'polarishub-
                        extended/tdcommons-cyp2d6-substrate-carbonmangels',
                        'polarishub-extended/tdcommons-cyp2d6-veith',
                        'polarishub-extended/tdcommons-cyp3a4-substrate-
                        carbonmangels', 'polarishub-extended/tdcommons-
                        cyp3a4-veith', 'polarishub-extended/tdcommons-dili',
                        'polarishub-extended/tdcommons-half-life-obach',
                        'polarishub-extended/tdcommons-herg', 'polarishub-
                        extended/tdcommons-hia-hou', 'polarishub-
                        extended/tdcommons-ld50-zhu', 'polarishub-
                        extended/tdcommons-lipophilicity-astrazeneca',
                        'polarishub-extended/tdcommons-pgp-broccatelli',
                        'polarishub-extended/tdcommons-ppbr-az', 'polarishub-
                        extended/tdcommons-solubility-aqsoldb', 'polarishub-
                        extended/tdcommons-vdss-lombardo', 'polarishub-
                        extended/vishrut64-cho-dna-expression-prediction-
                        dataset-task', 'polarishub-extended/vishrut64-rna-
                        expression-prediction-dataset-task', 'proteingym-
                        dms/A0A140D2T1_ZIKV_Sourisseau_2019', 'proteingym-
                        dms/A0A192B1T2_9HIV1_Haddox_2018', 'proteingym-
                        dms/A0A1I9GEU1_NEIME_Kennouche_2019', 'proteingym-
                        dms/A0A247D711_LISMN_Stadelmann_2021', 'proteingym-
                        dms/A0A2Z5U3Z0_9INFA_Doud_2016', 'proteingym-
                        dms/A0A2Z5U3Z0_9INFA_Wu_2014', 'proteingym-
                        dms/A4D664_9INFA_Soh_2019', 'proteingym-
                        dms/A4GRB6_PSEAI_Chen_2020', 'proteingym-
                        dms/A4_HUMAN_Seuma_2022', 'proteingym-
                        dms/A4_HUMAN_Seuma_2022_indels', 'proteingym-
                        dms/AACC1_PSEAI_Dandage_2018', 'proteingym-
                        dms/ACE2_HUMAN_Chan_2020', 'proteingym-
                        dms/ADRB2_HUMAN_Jones_2020', 'proteingym-
                        dms/AICDA_HUMAN_Gajula_2014_3cycles', 'proteingym-
                        dms/AMFR_HUMAN_Tsuboyama_2023_4G3O', 'proteingym-
                        dms/AMFR_HUMAN_Tsuboyama_2023_4G3O_indels',
                        'proteingym-dms/AMIE_PSEAE_Wrenbeck_2017',
                        'proteingym-dms/ANCSZ_Hobbs_2022', 'proteingym-
                        dms/ARGR_ECOLI_Tsuboyama_2023_1AOY', 'proteingym-
                        dms/ARGR_ECOLI_Tsuboyama_2023_1AOY_indels',
                        'proteingym-dms/B1LPA6_ECOSM_Russ_2020_indels',
                        'proteingym-dms/B2L11_HUMAN_Dutta_2010_binding-Mcl-1',
                        'proteingym-dms/BBC1_YEAST_Tsuboyama_2023_1TG0',
                        'proteingym-
                        dms/BBC1_YEAST_Tsuboyama_2023_1TG0_indels',
                        'proteingym-dms/BCHB_CHLTE_Tsuboyama_2023_2KRU',
                        'proteingym-
                        dms/BCHB_CHLTE_Tsuboyama_2023_2KRU_indels',
                        'proteingym-dms/BLAT_ECOLX_Deng_2012', 'proteingym-
                        dms/BLAT_ECOLX_Firnberg_2014', 'proteingym-
                        dms/BLAT_ECOLX_Gonzalez_2019_indels', 'proteingym-
                        dms/BLAT_ECOLX_Jacquier_2013', 'proteingym-
                        dms/BLAT_ECOLX_Stiffler_2015', 'proteingym-
                        dms/BRCA1_HUMAN_Findlay_2018', 'proteingym-
                        dms/BRCA2_HUMAN_Erwood_2022_HEK293T', 'proteingym-
                        dms/C6KNH7_9INFA_Lee_2018', 'proteingym-
                        dms/CALM1_HUMAN_Weile_2017', 'proteingym-
                        dms/CAPSD_AAV2S_Sinai_2021', 'proteingym-
                        dms/CAPSD_AAV2S_Sinai_2021_designed_indels',
                        'proteingym-
                        dms/CAPSD_AAV2S_Sinai_2021_library_indels',
                        'proteingym-dms/CAR11_HUMAN_Meitlis_2020_gof',
                        'proteingym-dms/CAR11_HUMAN_Meitlis_2020_lof',
                        'proteingym-dms/CAS9_STRP1_Spencer_2017_positive',
                        'proteingym-dms/CASP3_HUMAN_Roychowdhury_2020',
                        'proteingym-dms/CASP7_HUMAN_Roychowdhury_2020',
                        'proteingym-dms/CATR_CHLRE_Tsuboyama_2023_2AMI',
                        'proteingym-
                        dms/CATR_CHLRE_Tsuboyama_2023_2AMI_indels',
                        'proteingym-dms/CBPA2_HUMAN_Tsuboyama_2023_1O6X',
                        'proteingym-
                        dms/CBPA2_HUMAN_Tsuboyama_2023_1O6X_indels',
                        'proteingym-dms/CBS_HUMAN_Sun_2020', 'proteingym-
                        dms/CBX4_HUMAN_Tsuboyama_2023_2K28', 'proteingym-
                        dms/CBX4_HUMAN_Tsuboyama_2023_2K28_indels',
                        'proteingym-dms/CCDB_ECOLI_Adkar_2012', 'proteingym-
                        dms/CCDB_ECOLI_Tripathi_2016', 'proteingym-
                        dms/CCR5_HUMAN_Gill_2023', 'proteingym-
                        dms/CD19_HUMAN_Klesmith_2019_FMC_singles',
                        'proteingym-dms/CP2C9_HUMAN_Amorosi_2021_abundance',
                        'proteingym-dms/CP2C9_HUMAN_Amorosi_2021_activity',
                        'proteingym-dms/CSN4_MOUSE_Tsuboyama_2023_1UFM',
                        'proteingym-
                        dms/CSN4_MOUSE_Tsuboyama_2023_1UFM_indels',
                        'proteingym-dms/CUE1_YEAST_Tsuboyama_2023_2MYX',
                        'proteingym-
                        dms/CUE1_YEAST_Tsuboyama_2023_2MYX_indels',
                        'proteingym-dms/D7PM05_CLYGR_Somermeyer_2022',
                        'proteingym-dms/DLG4_HUMAN_Faure_2021', 'proteingym-
                        dms/DLG4_RAT_McLaughlin_2012', 'proteingym-
                        dms/DN7A_SACS2_Tsuboyama_2023_1JIC', 'proteingym-
                        dms/DN7A_SACS2_Tsuboyama_2023_1JIC_indels',
                        'proteingym-dms/DNJA1_HUMAN_Tsuboyama_2023_2LO1',
                        'proteingym-
                        dms/DNJA1_HUMAN_Tsuboyama_2023_2LO1_indels',
                        'proteingym-dms/DOCK1_MOUSE_Tsuboyama_2023_2M0Y',
                        'proteingym-
                        dms/DOCK1_MOUSE_Tsuboyama_2023_2M0Y_indels',
                        'proteingym-dms/DYR_ECOLI_Nguyen_2023', 'proteingym-
                        dms/DYR_ECOLI_Thompson_2019', 'proteingym-
                        dms/ENVZ_ECOLI_Ghose_2023', 'proteingym-
                        dms/ENV_HV1B9_DuenasDecamp_2016', 'proteingym-
                        dms/ENV_HV1BR_Haddox_2016', 'proteingym-
                        dms/EPHB2_HUMAN_Tsuboyama_2023_1F0M', 'proteingym-
                        dms/EPHB2_HUMAN_Tsuboyama_2023_1F0M_indels',
                        'proteingym-dms/ERBB2_HUMAN_Elazar_2016', 'proteingym-
                        dms/ESTA_BACSU_Nutschel_2020', 'proteingym-
                        dms/F7YBW8_MESOW_Aakre_2015', 'proteingym-
                        dms/F7YBW8_MESOW_Ding_2023', 'proteingym-
                        dms/FECA_ECOLI_Tsuboyama_2023_2D1U', 'proteingym-
                        dms/FECA_ECOLI_Tsuboyama_2023_2D1U_indels',
                        'proteingym-dms/FKBP3_HUMAN_Tsuboyama_2023_2KFV',
                        'proteingym-dms/GAL4_YEAST_Kitzman_2015', 'proteingym-
                        dms/GCN4_YEAST_Staller_2018', 'proteingym-
                        dms/GDIA_HUMAN_Silverstein_2021', 'proteingym-
                        dms/GFP_AEQVI_Sarkisyan_2016', 'proteingym-
                        dms/GLPA_HUMAN_Elazar_2016', 'proteingym-
                        dms/GRB2_HUMAN_Faure_2021', 'proteingym-
                        dms/HCP_LAMBD_Tsuboyama_2023_2L6Q', 'proteingym-
                        dms/HCP_LAMBD_Tsuboyama_2023_2L6Q_indels',
                        'proteingym-dms/HECD1_HUMAN_Tsuboyama_2023_3DKM',
                        'proteingym-
                        dms/HECD1_HUMAN_Tsuboyama_2023_3DKM_indels',
                        'proteingym-dms/HEM3_HUMAN_Loggerenberg_2023',
                        'proteingym-dms/HIS7_YEAST_Pokusaeva_2019',
                        'proteingym-dms/HIS7_YEAST_Pokusaeva_2019_indels',
                        'proteingym-dms/HMDH_HUMAN_Jiang_2019', 'proteingym-
                        dms/HSP82_YEAST_Cote-Hammarlof_2020_growth-H2O2',
                        'proteingym-dms/HSP82_YEAST_Flynn_2019', 'proteingym-
                        dms/HSP82_YEAST_Mishra_2016', 'proteingym-
                        dms/HXK4_HUMAN_Gersing_2022_activity', 'proteingym-
                        dms/HXK4_HUMAN_Gersing_2023_abundance', 'proteingym-
                        dms/I6TAH8_I68A0_Doud_2015', 'proteingym-
                        dms/IF1_ECOLI_Kelsic_2016', 'proteingym-
                        dms/ILF3_HUMAN_Tsuboyama_2023_2L33', 'proteingym-
                        dms/ILF3_HUMAN_Tsuboyama_2023_2L33_indels',
                        'proteingym-dms/ISDH_STAAW_Tsuboyama_2023_2LHR',
                        'proteingym-dms/KCNE1_HUMAN_Muhammad_2023_expression',
                        'proteingym-dms/KCNE1_HUMAN_Muhammad_2023_function',
                        'proteingym-dms/KCNH2_HUMAN_Kozek_2020', 'proteingym-
                        dms/KCNJ2_MOUSE_Coyote-Maestas_2022_function',
                        'proteingym-dms/KCNJ2_MOUSE_Coyote-
                        Maestas_2022_surface', 'proteingym-
                        dms/KCNJ2_MOUSE_Macdonald_2022_indels', 'proteingym-
                        dms/KKA2_KLEPN_Melnikov_2014', 'proteingym-
                        dms/LGK_LIPST_Klesmith_2015', 'proteingym-
                        dms/LYAM1_HUMAN_Elazar_2016', 'proteingym-
                        dms/MAFG_MOUSE_Tsuboyama_2023_1K1V', 'proteingym-
                        dms/MAFG_MOUSE_Tsuboyama_2023_1K1V_indels',
                        'proteingym-dms/MBD11_ARATH_Tsuboyama_2023_6ACV',
                        'proteingym-
                        dms/MBD11_ARATH_Tsuboyama_2023_6ACV_indels',
                        'proteingym-dms/MET_HUMAN_Estevam_2023', 'proteingym-
                        dms/MK01_HUMAN_Brenan_2016', 'proteingym-
                        dms/MLAC_ECOLI_MacRae_2023', 'proteingym-
                        dms/MSH2_HUMAN_Jia_2020', 'proteingym-
                        dms/MTH3_HAEAE_RockahShmuel_2015', 'proteingym-
                        dms/MTHR_HUMAN_Weile_2021', 'proteingym-
                        dms/MYO3_YEAST_Tsuboyama_2023_2BTT', 'proteingym-
                        dms/MYO3_YEAST_Tsuboyama_2023_2BTT_indels',
                        'proteingym-dms/NCAP_I34A1_Doud_2015', 'proteingym-
                        dms/NKX31_HUMAN_Tsuboyama_2023_2L9R', 'proteingym-
                        dms/NKX31_HUMAN_Tsuboyama_2023_2L9R_indels',
                        'proteingym-dms/NPC1_HUMAN_Erwood_2022_HEK293T',
                        'proteingym-dms/NPC1_HUMAN_Erwood_2022_RPE1',
                        'proteingym-dms/NRAM_I33A0_Jiang_2016', 'proteingym-
                        dms/NUD15_HUMAN_Suiter_2020', 'proteingym-
                        dms/NUSA_ECOLI_Tsuboyama_2023_1WCL', 'proteingym-
                        dms/NUSA_ECOLI_Tsuboyama_2023_1WCL_indels',
                        'proteingym-dms/NUSG_MYCTU_Tsuboyama_2023_2MI6',
                        'proteingym-
                        dms/NUSG_MYCTU_Tsuboyama_2023_2MI6_indels',
                        'proteingym-dms/OBSCN_HUMAN_Tsuboyama_2023_1V1C',
                        'proteingym-
                        dms/OBSCN_HUMAN_Tsuboyama_2023_1V1C_indels',
                        'proteingym-dms/ODP2_GEOSE_Tsuboyama_2023_1W4G',
                        'proteingym-
                        dms/ODP2_GEOSE_Tsuboyama_2023_1W4G_indels',
                        'proteingym-dms/OPSD_HUMAN_Wan_2019', 'proteingym-
                        dms/OTC_HUMAN_Lo_2023', 'proteingym-
                        dms/OTU7A_HUMAN_Tsuboyama_2023_2L2D', 'proteingym-
                        dms/OTU7A_HUMAN_Tsuboyama_2023_2L2D_indels',
                        'proteingym-dms/OXDA_RHOTO_Vanella_2023_activity',
                        'proteingym-dms/OXDA_RHOTO_Vanella_2023_expression',
                        'proteingym-
                        dms/P53_HUMAN_Giacomelli_2018_Null_Etoposide',
                        'proteingym-
                        dms/P53_HUMAN_Giacomelli_2018_Null_Nutlin',
                        'proteingym-dms/P53_HUMAN_Giacomelli_2018_WT_Nutlin',
                        'proteingym-dms/P53_HUMAN_Kotler_2018', 'proteingym-
                        dms/P53_HUMAN_Kotler_2018_indels', 'proteingym-
                        dms/P84126_THETH_Chan_2017', 'proteingym-
                        dms/PABP_YEAST_Melamed_2013', 'proteingym-
                        dms/PAI1_HUMAN_Huttinger_2021', 'proteingym-
                        dms/PA_I34A1_Wu_2015', 'proteingym-
                        dms/PHOT_CHLRE_Chen_2023', 'proteingym-
                        dms/PIN1_HUMAN_Tsuboyama_2023_1I6C', 'proteingym-
                        dms/PIN1_HUMAN_Tsuboyama_2023_1I6C_indels',
                        'proteingym-dms/PITX2_HUMAN_Tsuboyama_2023_2L7M',
                        'proteingym-
                        dms/PITX2_HUMAN_Tsuboyama_2023_2L7M_indels',
                        'proteingym-dms/PKN1_HUMAN_Tsuboyama_2023_1URF',
                        'proteingym-
                        dms/PKN1_HUMAN_Tsuboyama_2023_1URF_indels',
                        'proteingym-dms/POLG_CXB3N_Mattenberger_2021',
                        'proteingym-dms/POLG_DEN26_Suphatrakul_2023',
                        'proteingym-dms/POLG_HCVJF_Qi_2014', 'proteingym-
                        dms/POLG_PESV_Tsuboyama_2023_2MXD', 'proteingym-
                        dms/POLG_PESV_Tsuboyama_2023_2MXD_indels',
                        'proteingym-dms/PPARG_HUMAN_Majithia_2016',
                        'proteingym-dms/PPM1D_HUMAN_Miller_2022', 'proteingym-
                        dms/PR40A_HUMAN_Tsuboyama_2023_1UZC', 'proteingym-
                        dms/PR40A_HUMAN_Tsuboyama_2023_1UZC_indels',
                        'proteingym-dms/PRKN_HUMAN_Clausen_2023', 'proteingym-
                        dms/PSAE_PICP2_Tsuboyama_2023_1PSE', 'proteingym-
                        dms/PSAE_PICP2_Tsuboyama_2023_1PSE_indels',
                        'proteingym-dms/PTEN_HUMAN_Matreyek_2021',
                        'proteingym-dms/PTEN_HUMAN_Mighell_2018', 'proteingym-
                        dms/PTEN_HUMAN_Mighell_2018_indels', 'proteingym-
                        dms/Q2N0S5_9HIV1_Haddox_2018', 'proteingym-
                        dms/Q53Z42_HUMAN_McShan_2019_binding-TAPBPR',
                        'proteingym-dms/Q53Z42_HUMAN_McShan_2019_expression',
                        'proteingym-dms/Q59976_STRSQ_Romero_2015',
                        'proteingym-dms/Q6WV12_9MAXI_Somermeyer_2022',
                        'proteingym-dms/Q837P4_ENTFA_Meier_2023', 'proteingym-
                        dms/Q837P5_ENTFA_Meier_2023', 'proteingym-
                        dms/Q8EG35_SHEON_Campbell_2022_indels', 'proteingym-
                        dms/Q8WTC7_9CNID_Somermeyer_2022', 'proteingym-
                        dms/R1AB_SARS2_Flynn_2022', 'proteingym-
                        dms/RAD_ANTMA_Tsuboyama_2023_2CJJ', 'proteingym-
                        dms/RAD_ANTMA_Tsuboyama_2023_2CJJ_indels',
                        'proteingym-dms/RAF1_HUMAN_Zinkus-Boltz_2019',
                        'proteingym-dms/RASH_HUMAN_Bandaru_2017', 'proteingym-
                        dms/RASK_HUMAN_Weng_2022_abundance', 'proteingym-
                        dms/RASK_HUMAN_Weng_2022_binding-DARPin_K55',
                        'proteingym-dms/RBP1_HUMAN_Tsuboyama_2023_2KWH',
                        'proteingym-dms/RCD1_ARATH_Tsuboyama_2023_5OAO',
                        'proteingym-
                        dms/RCD1_ARATH_Tsuboyama_2023_5OAO_indels',
                        'proteingym-dms/RCRO_LAMBD_Tsuboyama_2023_1ORC',
                        'proteingym-dms/RD23A_HUMAN_Tsuboyama_2023_1IFY',
                        'proteingym-
                        dms/RD23A_HUMAN_Tsuboyama_2023_1IFY_indels',
                        'proteingym-dms/RDRP_I33A0_Li_2023', 'proteingym-
                        dms/REV_HV1H2_Fernandes_2016', 'proteingym-
                        dms/RFAH_ECOLI_Tsuboyama_2023_2LCL', 'proteingym-
                        dms/RL20_AQUAE_Tsuboyama_2023_1GYZ', 'proteingym-
                        dms/RL40A_YEAST_Mavor_2016', 'proteingym-
                        dms/RL40A_YEAST_Roscoe_2013', 'proteingym-
                        dms/RL40A_YEAST_Roscoe_2014', 'proteingym-
                        dms/RNC_ECOLI_Weeks_2023', 'proteingym-
                        dms/RPC1_BP434_Tsuboyama_2023_1R69', 'proteingym-
                        dms/RPC1_BP434_Tsuboyama_2023_1R69_indels',
                        'proteingym-dms/RPC1_LAMBD_Li_2019_high-expression',
                        'proteingym-dms/RPC1_LAMBD_Li_2019_low-expression',
                        'proteingym-dms/RS15_GEOSE_Tsuboyama_2023_1A32',
                        'proteingym-
                        dms/RS15_GEOSE_Tsuboyama_2023_1A32_indels',
                        'proteingym-dms/S22A1_HUMAN_Yee_2023_abundance',
                        'proteingym-
                        dms/S22A1_HUMAN_Yee_2023_abundance_indels',
                        'proteingym-dms/S22A1_HUMAN_Yee_2023_activity',
                        'proteingym-dms/S22A1_HUMAN_Yee_2023_activity_indels',
                        'proteingym-dms/SAV1_MOUSE_Tsuboyama_2023_2YSB',
                        'proteingym-
                        dms/SAV1_MOUSE_Tsuboyama_2023_2YSB_indels',
                        'proteingym-dms/SBI_STAAM_Tsuboyama_2023_2JVG',
                        'proteingym-dms/SC6A4_HUMAN_Young_2021', 'proteingym-
                        dms/SCIN_STAAR_Tsuboyama_2023_2QFF', 'proteingym-
                        dms/SCN5A_HUMAN_Glazer_2019', 'proteingym-
                        dms/SDA_BACSU_Tsuboyama_2023_1PV0', 'proteingym-
                        dms/SDA_BACSU_Tsuboyama_2023_1PV0_indels',
                        'proteingym-dms/SERC_HUMAN_Xie_2023', 'proteingym-
                        dms/SHOC2_HUMAN_Kwon_2022', 'proteingym-
                        dms/SOX30_HUMAN_Tsuboyama_2023_7JJK', 'proteingym-
                        dms/SOX30_HUMAN_Tsuboyama_2023_7JJK_indels',
                        'proteingym-dms/SPA_STAAU_Tsuboyama_2023_1LP1',
                        'proteingym-dms/SPG1_STRSG_Olson_2014', 'proteingym-
                        dms/SPG1_STRSG_Wu_2016', 'proteingym-
                        dms/SPG2_STRSG_Tsuboyama_2023_5UBS', 'proteingym-
                        dms/SPG2_STRSG_Tsuboyama_2023_5UBS_indels',
                        'proteingym-dms/SPIKE_SARS2_Starr_2020_binding',
                        'proteingym-dms/SPIKE_SARS2_Starr_2020_expression',
                        'proteingym-dms/SPTN1_CHICK_Tsuboyama_2023_1TUD',
                        'proteingym-
                        dms/SPTN1_CHICK_Tsuboyama_2023_1TUD_indels',
                        'proteingym-dms/SQSTM_MOUSE_Tsuboyama_2023_2RRU',
                        'proteingym-
                        dms/SQSTM_MOUSE_Tsuboyama_2023_2RRU_indels',
                        'proteingym-dms/SR43C_ARATH_Tsuboyama_2023_2N88',
                        'proteingym-
                        dms/SR43C_ARATH_Tsuboyama_2023_2N88_indels',
                        'proteingym-dms/SRBS1_HUMAN_Tsuboyama_2023_2O2W',
                        'proteingym-
                        dms/SRBS1_HUMAN_Tsuboyama_2023_2O2W_indels',
                        'proteingym-dms/SRC_HUMAN_Ahler_2019', 'proteingym-
                        dms/SRC_HUMAN_Chakraborty_2023_binding-DAS_25uM',
                        'proteingym-dms/SRC_HUMAN_Nguyen_2022', 'proteingym-
                        dms/SUMO1_HUMAN_Weile_2017', 'proteingym-
                        dms/SYUA_HUMAN_Newberry_2020', 'proteingym-
                        dms/TADBP_HUMAN_Bolognesi_2019', 'proteingym-
                        dms/TAT_HV1BR_Fernandes_2016', 'proteingym-
                        dms/TCRG1_MOUSE_Tsuboyama_2023_1E0L', 'proteingym-
                        dms/TCRG1_MOUSE_Tsuboyama_2023_1E0L_indels',
                        'proteingym-dms/THO1_YEAST_Tsuboyama_2023_2WQG',
                        'proteingym-
                        dms/THO1_YEAST_Tsuboyama_2023_2WQG_indels',
                        'proteingym-dms/TNKS2_HUMAN_Tsuboyama_2023_5JRT',
                        'proteingym-
                        dms/TNKS2_HUMAN_Tsuboyama_2023_5JRT_indels',
                        'proteingym-dms/TPK1_HUMAN_Weile_2017', 'proteingym-
                        dms/TPMT_HUMAN_Matreyek_2018', 'proteingym-
                        dms/TPOR_HUMAN_Bridgford_2020', 'proteingym-
                        dms/TRPC_SACS2_Chan_2017', 'proteingym-
                        dms/TRPC_THEMA_Chan_2017', 'proteingym-
                        dms/UBC9_HUMAN_Weile_2017', 'proteingym-
                        dms/UBE4B_HUMAN_Tsuboyama_2023_3L1X', 'proteingym-
                        dms/UBE4B_HUMAN_Tsuboyama_2023_3L1X_indels',
                        'proteingym-dms/UBE4B_MOUSE_Starita_2013',
                        'proteingym-dms/UBR5_HUMAN_Tsuboyama_2023_1I2T',
                        'proteingym-
                        dms/UBR5_HUMAN_Tsuboyama_2023_1I2T_indels',
                        'proteingym-dms/VG08_BPP22_Tsuboyama_2023_2GP8',
                        'proteingym-
                        dms/VG08_BPP22_Tsuboyama_2023_2GP8_indels',
                        'proteingym-dms/VILI_CHICK_Tsuboyama_2023_1YU5',
                        'proteingym-
                        dms/VILI_CHICK_Tsuboyama_2023_1YU5_indels',
                        'proteingym-dms/VKOR1_HUMAN_Chiasson_2020_abundance',
                        'proteingym-dms/VKOR1_HUMAN_Chiasson_2020_activity',
                        'proteingym-dms/VRPI_BPT7_Tsuboyama_2023_2WNM',
                        'proteingym-dms/VRPI_BPT7_Tsuboyama_2023_2WNM_indels',
                        'proteingym-dms/YAIA_ECOLI_Tsuboyama_2023_2KVT',
                        'proteingym-dms/YAP1_HUMAN_Araya_2012', 'proteingym-
                        dms/YNZC_BACSU_Tsuboyama_2023_2JVD', 'proteingym-
                        dms/YNZC_BACSU_Tsuboyama_2023_2JVD_indels']
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
                        options: ['kaggle/histopathologic-cancer-detection',
                        'kaggle/hms-harmful-brain-activity-classification',
                        'kaggle/osic-pulmonary-fibrosis-progression',
                        'kaggle/ranzcr-clip-catheter-line-classification',
                        'kaggle/rsna-miccai-brain-tumor-radiogenomic-
                        classification', 'kaggle/stanford-covid-vaccine',
                        'kaggle/uw-madison-gi-tract-image-segmentation',
                        'kaggle/ventilator-pressure-prediction',
                        'kaggle_big/bms-molecular-translation',
                        'kaggle_big/rsna-2022-cervical-spine-fracture-
                        detection', 'kaggle_big/rsna-breast-cancer-detection',
                        'kaggle_big/siim-covid19-detection', 'kaggle_big/siim-
                        isic-melanoma-classification', 'kaggle_big/vinbigdata-
                        chest-xray-abnormalities-detection', 'manual/open-
                        problems-cell-cell-communication-ligand-target',
                        'manual/open-problems-label-projection', 'manual/open-
                        problems-predict-modality', 'manual/open-problems-
                        single-cell-perturbations', 'manual/open-problems-
                        spatially-variable-genes', 'polarishub/polaris-adme-
                        fang-hclint-1', 'polarishub/polaris-adme-fang-hppb-1',
                        'polarishub/polaris-adme-fang-perm-1',
                        'polarishub/polaris-adme-fang-r-1',
                        'polarishub/polaris-adme-fang-rclint-1',
                        'polarishub/polaris-adme-fang-rppb-1',
                        'polarishub/polaris-adme-fang-solu-1',
                        'polarishub/polaris-pkis1-egfr-wt-mut-c-1',
                        'polarishub/polaris-pkis1-egfr-wt-mut-r-1',
                        'polarishub/polaris-pkis1-kit-wt-mut-c-1',
                        'polarishub/polaris-pkis1-kit-wt-mut-r-1',
                        'polarishub/polaris-pkis1-ret-wt-mut-c-1',
                        'polarishub/polaris-pkis1-ret-wt-mut-r-1',
                        'polarishub/polaris-pkis2-egfr-wt-c-1',
                        'polarishub/polaris-pkis2-egfr-wt-r-1',
                        'polarishub/polaris-pkis2-kit-wt-c-1',
                        'polarishub/polaris-pkis2-kit-wt-r-1',
                        'polarishub/polaris-pkis2-lok-slk-c-1',
                        'polarishub/polaris-pkis2-lok-slk-r-1',
                        'polarishub/polaris-pkis2-ret-wt-c-1',
                        'polarishub/polaris-pkis2-ret-wt-r-1',
                        'polarishub/tdcommons-bbb-martins',
                        'polarishub/tdcommons-bioavailability-ma',
                        'polarishub/tdcommons-caco2-wang',
                        'polarishub/tdcommons-clearance-hepatocyte-az',
                        'polarishub/tdcommons-cyp2d6-substrate-carbonmangels',
                        'polarishub/tdcommons-herg', 'polarishub/tdcommons-
                        hia-hou', 'polarishub/tdcommons-lipophilicity-
                        astrazeneca', 'polarishub/tdcommons-pgp-broccatelli',
                        'polarishub-extended/adaptyv-bio-egfr-binders-binary-
                        cls-v0', 'polarishub-extended/biogen-adme-fang-hclint-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-hppb-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-perm-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-rclint-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-rppb-
                        reg-v1', 'polarishub-extended/biogen-adme-fang-solu-
                        reg-v1', 'polarishub-extended/graphium-l1000-mcf7-v1',
                        'polarishub-extended/graphium-l1000-vcap-v1',
                        'polarishub-extended/graphium-pcba-1328-1564k-v1',
                        'polarishub-extended/graphium-qm9-v1', 'polarishub-
                        extended/graphium-tox21-v1', 'polarishub-
                        extended/graphium-zinc12k-v1', 'polarishub-
                        extended/novartis-adme-novartis-cyp3a4-cls',
                        'polarishub-extended/novartis-adme-novartis-
                        cyp3a4-reg', 'polarishub-extended/polaris-adme-fang-
                        hclint-1', 'polarishub-extended/polaris-adme-fang-
                        hppb-1', 'polarishub-extended/polaris-adme-fang-
                        perm-1', 'polarishub-extended/polaris-adme-fang-r-1',
                        'polarishub-extended/polaris-adme-fang-rclint-1',
                        'polarishub-extended/polaris-adme-fang-rppb-1',
                        'polarishub-extended/polaris-adme-fang-solu-1',
                        'polarishub-extended/polaris-hello-world-benchmark',
                        'polarishub-extended/polaris-
                        molprop-250k-leadlike-r-1', 'polarishub-
                        extended/polaris-molprop-250k-r-1', 'polarishub-
                        extended/polaris-molprop-250k-reg-v2', 'polarishub-
                        extended/polaris-molprop-leadlike-250k-reg-v2',
                        'polarishub-extended/polaris-pkis1-egfr-wt-mut-c-1',
                        'polarishub-extended/polaris-pkis1-egfr-wt-mut-r-1',
                        'polarishub-extended/polaris-pkis1-kit-wt-mut-c-1',
                        'polarishub-extended/polaris-pkis1-kit-wt-mut-r-1',
                        'polarishub-extended/polaris-pkis1-ret-wt-mut-c-1',
                        'polarishub-extended/polaris-pkis1-ret-wt-mut-r-1',
                        'polarishub-extended/polaris-pkis2-egfr-wt-c-1',
                        'polarishub-extended/polaris-pkis2-egfr-wt-r-1',
                        'polarishub-extended/polaris-pkis2-egfr-wt-reg-v2',
                        'polarishub-extended/polaris-pkis2-kit-wt-c-1',
                        'polarishub-extended/polaris-pkis2-kit-wt-cls-v2',
                        'polarishub-extended/polaris-pkis2-kit-wt-r-1',
                        'polarishub-extended/polaris-pkis2-kit-wt-reg-v2',
                        'polarishub-extended/polaris-pkis2-lok-slk-c-1',
                        'polarishub-extended/polaris-pkis2-lok-slk-cls-v2',
                        'polarishub-extended/polaris-pkis2-lok-slk-r-1',
                        'polarishub-extended/polaris-pkis2-ret-wt-c-1',
                        'polarishub-extended/polaris-pkis2-ret-wt-cls-v2',
                        'polarishub-extended/polaris-pkis2-ret-wt-r-1',
                        'polarishub-extended/polaris-pkis2-ret-wt-reg-v2',
                        'polarishub-extended/polaris-posebusters-v1',
                        'polarishub-extended/recursion-rxrx-compound-gene-
                        activity-benchmark', 'polarishub-extended/tdcommons-
                        ames', 'polarishub-extended/tdcommons-bbb-martins',
                        'polarishub-extended/tdcommons-bioavailability-ma',
                        'polarishub-extended/tdcommons-caco2-wang',
                        'polarishub-extended/tdcommons-clearance-hepatocyte-
                        az', 'polarishub-extended/tdcommons-clearance-
                        microsome-az', 'polarishub-extended/tdcommons-
                        cyp2c9-substrate-carbonmangels', 'polarishub-
                        extended/tdcommons-cyp2c9-veith', 'polarishub-
                        extended/tdcommons-cyp2d6-substrate-carbonmangels',
                        'polarishub-extended/tdcommons-cyp2d6-veith',
                        'polarishub-extended/tdcommons-cyp3a4-substrate-
                        carbonmangels', 'polarishub-extended/tdcommons-
                        cyp3a4-veith', 'polarishub-extended/tdcommons-dili',
                        'polarishub-extended/tdcommons-half-life-obach',
                        'polarishub-extended/tdcommons-herg', 'polarishub-
                        extended/tdcommons-hia-hou', 'polarishub-
                        extended/tdcommons-ld50-zhu', 'polarishub-
                        extended/tdcommons-lipophilicity-astrazeneca',
                        'polarishub-extended/tdcommons-pgp-broccatelli',
                        'polarishub-extended/tdcommons-ppbr-az', 'polarishub-
                        extended/tdcommons-solubility-aqsoldb', 'polarishub-
                        extended/tdcommons-vdss-lombardo', 'polarishub-
                        extended/vishrut64-cho-dna-expression-prediction-
                        dataset-task', 'polarishub-extended/vishrut64-rna-
                        expression-prediction-dataset-task', 'proteingym-
                        dms/A0A140D2T1_ZIKV_Sourisseau_2019', 'proteingym-
                        dms/A0A192B1T2_9HIV1_Haddox_2018', 'proteingym-
                        dms/A0A1I9GEU1_NEIME_Kennouche_2019', 'proteingym-
                        dms/A0A247D711_LISMN_Stadelmann_2021', 'proteingym-
                        dms/A0A2Z5U3Z0_9INFA_Doud_2016', 'proteingym-
                        dms/A0A2Z5U3Z0_9INFA_Wu_2014', 'proteingym-
                        dms/A4D664_9INFA_Soh_2019', 'proteingym-
                        dms/A4GRB6_PSEAI_Chen_2020', 'proteingym-
                        dms/A4_HUMAN_Seuma_2022', 'proteingym-
                        dms/A4_HUMAN_Seuma_2022_indels', 'proteingym-
                        dms/AACC1_PSEAI_Dandage_2018', 'proteingym-
                        dms/ACE2_HUMAN_Chan_2020', 'proteingym-
                        dms/ADRB2_HUMAN_Jones_2020', 'proteingym-
                        dms/AICDA_HUMAN_Gajula_2014_3cycles', 'proteingym-
                        dms/AMFR_HUMAN_Tsuboyama_2023_4G3O', 'proteingym-
                        dms/AMFR_HUMAN_Tsuboyama_2023_4G3O_indels',
                        'proteingym-dms/AMIE_PSEAE_Wrenbeck_2017',
                        'proteingym-dms/ANCSZ_Hobbs_2022', 'proteingym-
                        dms/ARGR_ECOLI_Tsuboyama_2023_1AOY', 'proteingym-
                        dms/ARGR_ECOLI_Tsuboyama_2023_1AOY_indels',
                        'proteingym-dms/B1LPA6_ECOSM_Russ_2020_indels',
                        'proteingym-dms/B2L11_HUMAN_Dutta_2010_binding-Mcl-1',
                        'proteingym-dms/BBC1_YEAST_Tsuboyama_2023_1TG0',
                        'proteingym-
                        dms/BBC1_YEAST_Tsuboyama_2023_1TG0_indels',
                        'proteingym-dms/BCHB_CHLTE_Tsuboyama_2023_2KRU',
                        'proteingym-
                        dms/BCHB_CHLTE_Tsuboyama_2023_2KRU_indels',
                        'proteingym-dms/BLAT_ECOLX_Deng_2012', 'proteingym-
                        dms/BLAT_ECOLX_Firnberg_2014', 'proteingym-
                        dms/BLAT_ECOLX_Gonzalez_2019_indels', 'proteingym-
                        dms/BLAT_ECOLX_Jacquier_2013', 'proteingym-
                        dms/BLAT_ECOLX_Stiffler_2015', 'proteingym-
                        dms/BRCA1_HUMAN_Findlay_2018', 'proteingym-
                        dms/BRCA2_HUMAN_Erwood_2022_HEK293T', 'proteingym-
                        dms/C6KNH7_9INFA_Lee_2018', 'proteingym-
                        dms/CALM1_HUMAN_Weile_2017', 'proteingym-
                        dms/CAPSD_AAV2S_Sinai_2021', 'proteingym-
                        dms/CAPSD_AAV2S_Sinai_2021_designed_indels',
                        'proteingym-
                        dms/CAPSD_AAV2S_Sinai_2021_library_indels',
                        'proteingym-dms/CAR11_HUMAN_Meitlis_2020_gof',
                        'proteingym-dms/CAR11_HUMAN_Meitlis_2020_lof',
                        'proteingym-dms/CAS9_STRP1_Spencer_2017_positive',
                        'proteingym-dms/CASP3_HUMAN_Roychowdhury_2020',
                        'proteingym-dms/CASP7_HUMAN_Roychowdhury_2020',
                        'proteingym-dms/CATR_CHLRE_Tsuboyama_2023_2AMI',
                        'proteingym-
                        dms/CATR_CHLRE_Tsuboyama_2023_2AMI_indels',
                        'proteingym-dms/CBPA2_HUMAN_Tsuboyama_2023_1O6X',
                        'proteingym-
                        dms/CBPA2_HUMAN_Tsuboyama_2023_1O6X_indels',
                        'proteingym-dms/CBS_HUMAN_Sun_2020', 'proteingym-
                        dms/CBX4_HUMAN_Tsuboyama_2023_2K28', 'proteingym-
                        dms/CBX4_HUMAN_Tsuboyama_2023_2K28_indels',
                        'proteingym-dms/CCDB_ECOLI_Adkar_2012', 'proteingym-
                        dms/CCDB_ECOLI_Tripathi_2016', 'proteingym-
                        dms/CCR5_HUMAN_Gill_2023', 'proteingym-
                        dms/CD19_HUMAN_Klesmith_2019_FMC_singles',
                        'proteingym-dms/CP2C9_HUMAN_Amorosi_2021_abundance',
                        'proteingym-dms/CP2C9_HUMAN_Amorosi_2021_activity',
                        'proteingym-dms/CSN4_MOUSE_Tsuboyama_2023_1UFM',
                        'proteingym-
                        dms/CSN4_MOUSE_Tsuboyama_2023_1UFM_indels',
                        'proteingym-dms/CUE1_YEAST_Tsuboyama_2023_2MYX',
                        'proteingym-
                        dms/CUE1_YEAST_Tsuboyama_2023_2MYX_indels',
                        'proteingym-dms/D7PM05_CLYGR_Somermeyer_2022',
                        'proteingym-dms/DLG4_HUMAN_Faure_2021', 'proteingym-
                        dms/DLG4_RAT_McLaughlin_2012', 'proteingym-
                        dms/DN7A_SACS2_Tsuboyama_2023_1JIC', 'proteingym-
                        dms/DN7A_SACS2_Tsuboyama_2023_1JIC_indels',
                        'proteingym-dms/DNJA1_HUMAN_Tsuboyama_2023_2LO1',
                        'proteingym-
                        dms/DNJA1_HUMAN_Tsuboyama_2023_2LO1_indels',
                        'proteingym-dms/DOCK1_MOUSE_Tsuboyama_2023_2M0Y',
                        'proteingym-
                        dms/DOCK1_MOUSE_Tsuboyama_2023_2M0Y_indels',
                        'proteingym-dms/DYR_ECOLI_Nguyen_2023', 'proteingym-
                        dms/DYR_ECOLI_Thompson_2019', 'proteingym-
                        dms/ENVZ_ECOLI_Ghose_2023', 'proteingym-
                        dms/ENV_HV1B9_DuenasDecamp_2016', 'proteingym-
                        dms/ENV_HV1BR_Haddox_2016', 'proteingym-
                        dms/EPHB2_HUMAN_Tsuboyama_2023_1F0M', 'proteingym-
                        dms/EPHB2_HUMAN_Tsuboyama_2023_1F0M_indels',
                        'proteingym-dms/ERBB2_HUMAN_Elazar_2016', 'proteingym-
                        dms/ESTA_BACSU_Nutschel_2020', 'proteingym-
                        dms/F7YBW8_MESOW_Aakre_2015', 'proteingym-
                        dms/F7YBW8_MESOW_Ding_2023', 'proteingym-
                        dms/FECA_ECOLI_Tsuboyama_2023_2D1U', 'proteingym-
                        dms/FECA_ECOLI_Tsuboyama_2023_2D1U_indels',
                        'proteingym-dms/FKBP3_HUMAN_Tsuboyama_2023_2KFV',
                        'proteingym-dms/GAL4_YEAST_Kitzman_2015', 'proteingym-
                        dms/GCN4_YEAST_Staller_2018', 'proteingym-
                        dms/GDIA_HUMAN_Silverstein_2021', 'proteingym-
                        dms/GFP_AEQVI_Sarkisyan_2016', 'proteingym-
                        dms/GLPA_HUMAN_Elazar_2016', 'proteingym-
                        dms/GRB2_HUMAN_Faure_2021', 'proteingym-
                        dms/HCP_LAMBD_Tsuboyama_2023_2L6Q', 'proteingym-
                        dms/HCP_LAMBD_Tsuboyama_2023_2L6Q_indels',
                        'proteingym-dms/HECD1_HUMAN_Tsuboyama_2023_3DKM',
                        'proteingym-
                        dms/HECD1_HUMAN_Tsuboyama_2023_3DKM_indels',
                        'proteingym-dms/HEM3_HUMAN_Loggerenberg_2023',
                        'proteingym-dms/HIS7_YEAST_Pokusaeva_2019',
                        'proteingym-dms/HIS7_YEAST_Pokusaeva_2019_indels',
                        'proteingym-dms/HMDH_HUMAN_Jiang_2019', 'proteingym-
                        dms/HSP82_YEAST_Cote-Hammarlof_2020_growth-H2O2',
                        'proteingym-dms/HSP82_YEAST_Flynn_2019', 'proteingym-
                        dms/HSP82_YEAST_Mishra_2016', 'proteingym-
                        dms/HXK4_HUMAN_Gersing_2022_activity', 'proteingym-
                        dms/HXK4_HUMAN_Gersing_2023_abundance', 'proteingym-
                        dms/I6TAH8_I68A0_Doud_2015', 'proteingym-
                        dms/IF1_ECOLI_Kelsic_2016', 'proteingym-
                        dms/ILF3_HUMAN_Tsuboyama_2023_2L33', 'proteingym-
                        dms/ILF3_HUMAN_Tsuboyama_2023_2L33_indels',
                        'proteingym-dms/ISDH_STAAW_Tsuboyama_2023_2LHR',
                        'proteingym-dms/KCNE1_HUMAN_Muhammad_2023_expression',
                        'proteingym-dms/KCNE1_HUMAN_Muhammad_2023_function',
                        'proteingym-dms/KCNH2_HUMAN_Kozek_2020', 'proteingym-
                        dms/KCNJ2_MOUSE_Coyote-Maestas_2022_function',
                        'proteingym-dms/KCNJ2_MOUSE_Coyote-
                        Maestas_2022_surface', 'proteingym-
                        dms/KCNJ2_MOUSE_Macdonald_2022_indels', 'proteingym-
                        dms/KKA2_KLEPN_Melnikov_2014', 'proteingym-
                        dms/LGK_LIPST_Klesmith_2015', 'proteingym-
                        dms/LYAM1_HUMAN_Elazar_2016', 'proteingym-
                        dms/MAFG_MOUSE_Tsuboyama_2023_1K1V', 'proteingym-
                        dms/MAFG_MOUSE_Tsuboyama_2023_1K1V_indels',
                        'proteingym-dms/MBD11_ARATH_Tsuboyama_2023_6ACV',
                        'proteingym-
                        dms/MBD11_ARATH_Tsuboyama_2023_6ACV_indels',
                        'proteingym-dms/MET_HUMAN_Estevam_2023', 'proteingym-
                        dms/MK01_HUMAN_Brenan_2016', 'proteingym-
                        dms/MLAC_ECOLI_MacRae_2023', 'proteingym-
                        dms/MSH2_HUMAN_Jia_2020', 'proteingym-
                        dms/MTH3_HAEAE_RockahShmuel_2015', 'proteingym-
                        dms/MTHR_HUMAN_Weile_2021', 'proteingym-
                        dms/MYO3_YEAST_Tsuboyama_2023_2BTT', 'proteingym-
                        dms/MYO3_YEAST_Tsuboyama_2023_2BTT_indels',
                        'proteingym-dms/NCAP_I34A1_Doud_2015', 'proteingym-
                        dms/NKX31_HUMAN_Tsuboyama_2023_2L9R', 'proteingym-
                        dms/NKX31_HUMAN_Tsuboyama_2023_2L9R_indels',
                        'proteingym-dms/NPC1_HUMAN_Erwood_2022_HEK293T',
                        'proteingym-dms/NPC1_HUMAN_Erwood_2022_RPE1',
                        'proteingym-dms/NRAM_I33A0_Jiang_2016', 'proteingym-
                        dms/NUD15_HUMAN_Suiter_2020', 'proteingym-
                        dms/NUSA_ECOLI_Tsuboyama_2023_1WCL', 'proteingym-
                        dms/NUSA_ECOLI_Tsuboyama_2023_1WCL_indels',
                        'proteingym-dms/NUSG_MYCTU_Tsuboyama_2023_2MI6',
                        'proteingym-
                        dms/NUSG_MYCTU_Tsuboyama_2023_2MI6_indels',
                        'proteingym-dms/OBSCN_HUMAN_Tsuboyama_2023_1V1C',
                        'proteingym-
                        dms/OBSCN_HUMAN_Tsuboyama_2023_1V1C_indels',
                        'proteingym-dms/ODP2_GEOSE_Tsuboyama_2023_1W4G',
                        'proteingym-
                        dms/ODP2_GEOSE_Tsuboyama_2023_1W4G_indels',
                        'proteingym-dms/OPSD_HUMAN_Wan_2019', 'proteingym-
                        dms/OTC_HUMAN_Lo_2023', 'proteingym-
                        dms/OTU7A_HUMAN_Tsuboyama_2023_2L2D', 'proteingym-
                        dms/OTU7A_HUMAN_Tsuboyama_2023_2L2D_indels',
                        'proteingym-dms/OXDA_RHOTO_Vanella_2023_activity',
                        'proteingym-dms/OXDA_RHOTO_Vanella_2023_expression',
                        'proteingym-
                        dms/P53_HUMAN_Giacomelli_2018_Null_Etoposide',
                        'proteingym-
                        dms/P53_HUMAN_Giacomelli_2018_Null_Nutlin',
                        'proteingym-dms/P53_HUMAN_Giacomelli_2018_WT_Nutlin',
                        'proteingym-dms/P53_HUMAN_Kotler_2018', 'proteingym-
                        dms/P53_HUMAN_Kotler_2018_indels', 'proteingym-
                        dms/P84126_THETH_Chan_2017', 'proteingym-
                        dms/PABP_YEAST_Melamed_2013', 'proteingym-
                        dms/PAI1_HUMAN_Huttinger_2021', 'proteingym-
                        dms/PA_I34A1_Wu_2015', 'proteingym-
                        dms/PHOT_CHLRE_Chen_2023', 'proteingym-
                        dms/PIN1_HUMAN_Tsuboyama_2023_1I6C', 'proteingym-
                        dms/PIN1_HUMAN_Tsuboyama_2023_1I6C_indels',
                        'proteingym-dms/PITX2_HUMAN_Tsuboyama_2023_2L7M',
                        'proteingym-
                        dms/PITX2_HUMAN_Tsuboyama_2023_2L7M_indels',
                        'proteingym-dms/PKN1_HUMAN_Tsuboyama_2023_1URF',
                        'proteingym-
                        dms/PKN1_HUMAN_Tsuboyama_2023_1URF_indels',
                        'proteingym-dms/POLG_CXB3N_Mattenberger_2021',
                        'proteingym-dms/POLG_DEN26_Suphatrakul_2023',
                        'proteingym-dms/POLG_HCVJF_Qi_2014', 'proteingym-
                        dms/POLG_PESV_Tsuboyama_2023_2MXD', 'proteingym-
                        dms/POLG_PESV_Tsuboyama_2023_2MXD_indels',
                        'proteingym-dms/PPARG_HUMAN_Majithia_2016',
                        'proteingym-dms/PPM1D_HUMAN_Miller_2022', 'proteingym-
                        dms/PR40A_HUMAN_Tsuboyama_2023_1UZC', 'proteingym-
                        dms/PR40A_HUMAN_Tsuboyama_2023_1UZC_indels',
                        'proteingym-dms/PRKN_HUMAN_Clausen_2023', 'proteingym-
                        dms/PSAE_PICP2_Tsuboyama_2023_1PSE', 'proteingym-
                        dms/PSAE_PICP2_Tsuboyama_2023_1PSE_indels',
                        'proteingym-dms/PTEN_HUMAN_Matreyek_2021',
                        'proteingym-dms/PTEN_HUMAN_Mighell_2018', 'proteingym-
                        dms/PTEN_HUMAN_Mighell_2018_indels', 'proteingym-
                        dms/Q2N0S5_9HIV1_Haddox_2018', 'proteingym-
                        dms/Q53Z42_HUMAN_McShan_2019_binding-TAPBPR',
                        'proteingym-dms/Q53Z42_HUMAN_McShan_2019_expression',
                        'proteingym-dms/Q59976_STRSQ_Romero_2015',
                        'proteingym-dms/Q6WV12_9MAXI_Somermeyer_2022',
                        'proteingym-dms/Q837P4_ENTFA_Meier_2023', 'proteingym-
                        dms/Q837P5_ENTFA_Meier_2023', 'proteingym-
                        dms/Q8EG35_SHEON_Campbell_2022_indels', 'proteingym-
                        dms/Q8WTC7_9CNID_Somermeyer_2022', 'proteingym-
                        dms/R1AB_SARS2_Flynn_2022', 'proteingym-
                        dms/RAD_ANTMA_Tsuboyama_2023_2CJJ', 'proteingym-
                        dms/RAD_ANTMA_Tsuboyama_2023_2CJJ_indels',
                        'proteingym-dms/RAF1_HUMAN_Zinkus-Boltz_2019',
                        'proteingym-dms/RASH_HUMAN_Bandaru_2017', 'proteingym-
                        dms/RASK_HUMAN_Weng_2022_abundance', 'proteingym-
                        dms/RASK_HUMAN_Weng_2022_binding-DARPin_K55',
                        'proteingym-dms/RBP1_HUMAN_Tsuboyama_2023_2KWH',
                        'proteingym-dms/RCD1_ARATH_Tsuboyama_2023_5OAO',
                        'proteingym-
                        dms/RCD1_ARATH_Tsuboyama_2023_5OAO_indels',
                        'proteingym-dms/RCRO_LAMBD_Tsuboyama_2023_1ORC',
                        'proteingym-dms/RD23A_HUMAN_Tsuboyama_2023_1IFY',
                        'proteingym-
                        dms/RD23A_HUMAN_Tsuboyama_2023_1IFY_indels',
                        'proteingym-dms/RDRP_I33A0_Li_2023', 'proteingym-
                        dms/REV_HV1H2_Fernandes_2016', 'proteingym-
                        dms/RFAH_ECOLI_Tsuboyama_2023_2LCL', 'proteingym-
                        dms/RL20_AQUAE_Tsuboyama_2023_1GYZ', 'proteingym-
                        dms/RL40A_YEAST_Mavor_2016', 'proteingym-
                        dms/RL40A_YEAST_Roscoe_2013', 'proteingym-
                        dms/RL40A_YEAST_Roscoe_2014', 'proteingym-
                        dms/RNC_ECOLI_Weeks_2023', 'proteingym-
                        dms/RPC1_BP434_Tsuboyama_2023_1R69', 'proteingym-
                        dms/RPC1_BP434_Tsuboyama_2023_1R69_indels',
                        'proteingym-dms/RPC1_LAMBD_Li_2019_high-expression',
                        'proteingym-dms/RPC1_LAMBD_Li_2019_low-expression',
                        'proteingym-dms/RS15_GEOSE_Tsuboyama_2023_1A32',
                        'proteingym-
                        dms/RS15_GEOSE_Tsuboyama_2023_1A32_indels',
                        'proteingym-dms/S22A1_HUMAN_Yee_2023_abundance',
                        'proteingym-
                        dms/S22A1_HUMAN_Yee_2023_abundance_indels',
                        'proteingym-dms/S22A1_HUMAN_Yee_2023_activity',
                        'proteingym-dms/S22A1_HUMAN_Yee_2023_activity_indels',
                        'proteingym-dms/SAV1_MOUSE_Tsuboyama_2023_2YSB',
                        'proteingym-
                        dms/SAV1_MOUSE_Tsuboyama_2023_2YSB_indels',
                        'proteingym-dms/SBI_STAAM_Tsuboyama_2023_2JVG',
                        'proteingym-dms/SC6A4_HUMAN_Young_2021', 'proteingym-
                        dms/SCIN_STAAR_Tsuboyama_2023_2QFF', 'proteingym-
                        dms/SCN5A_HUMAN_Glazer_2019', 'proteingym-
                        dms/SDA_BACSU_Tsuboyama_2023_1PV0', 'proteingym-
                        dms/SDA_BACSU_Tsuboyama_2023_1PV0_indels',
                        'proteingym-dms/SERC_HUMAN_Xie_2023', 'proteingym-
                        dms/SHOC2_HUMAN_Kwon_2022', 'proteingym-
                        dms/SOX30_HUMAN_Tsuboyama_2023_7JJK', 'proteingym-
                        dms/SOX30_HUMAN_Tsuboyama_2023_7JJK_indels',
                        'proteingym-dms/SPA_STAAU_Tsuboyama_2023_1LP1',
                        'proteingym-dms/SPG1_STRSG_Olson_2014', 'proteingym-
                        dms/SPG1_STRSG_Wu_2016', 'proteingym-
                        dms/SPG2_STRSG_Tsuboyama_2023_5UBS', 'proteingym-
                        dms/SPG2_STRSG_Tsuboyama_2023_5UBS_indels',
                        'proteingym-dms/SPIKE_SARS2_Starr_2020_binding',
                        'proteingym-dms/SPIKE_SARS2_Starr_2020_expression',
                        'proteingym-dms/SPTN1_CHICK_Tsuboyama_2023_1TUD',
                        'proteingym-
                        dms/SPTN1_CHICK_Tsuboyama_2023_1TUD_indels',
                        'proteingym-dms/SQSTM_MOUSE_Tsuboyama_2023_2RRU',
                        'proteingym-
                        dms/SQSTM_MOUSE_Tsuboyama_2023_2RRU_indels',
                        'proteingym-dms/SR43C_ARATH_Tsuboyama_2023_2N88',
                        'proteingym-
                        dms/SR43C_ARATH_Tsuboyama_2023_2N88_indels',
                        'proteingym-dms/SRBS1_HUMAN_Tsuboyama_2023_2O2W',
                        'proteingym-
                        dms/SRBS1_HUMAN_Tsuboyama_2023_2O2W_indels',
                        'proteingym-dms/SRC_HUMAN_Ahler_2019', 'proteingym-
                        dms/SRC_HUMAN_Chakraborty_2023_binding-DAS_25uM',
                        'proteingym-dms/SRC_HUMAN_Nguyen_2022', 'proteingym-
                        dms/SUMO1_HUMAN_Weile_2017', 'proteingym-
                        dms/SYUA_HUMAN_Newberry_2020', 'proteingym-
                        dms/TADBP_HUMAN_Bolognesi_2019', 'proteingym-
                        dms/TAT_HV1BR_Fernandes_2016', 'proteingym-
                        dms/TCRG1_MOUSE_Tsuboyama_2023_1E0L', 'proteingym-
                        dms/TCRG1_MOUSE_Tsuboyama_2023_1E0L_indels',
                        'proteingym-dms/THO1_YEAST_Tsuboyama_2023_2WQG',
                        'proteingym-
                        dms/THO1_YEAST_Tsuboyama_2023_2WQG_indels',
                        'proteingym-dms/TNKS2_HUMAN_Tsuboyama_2023_5JRT',
                        'proteingym-
                        dms/TNKS2_HUMAN_Tsuboyama_2023_5JRT_indels',
                        'proteingym-dms/TPK1_HUMAN_Weile_2017', 'proteingym-
                        dms/TPMT_HUMAN_Matreyek_2018', 'proteingym-
                        dms/TPOR_HUMAN_Bridgford_2020', 'proteingym-
                        dms/TRPC_SACS2_Chan_2017', 'proteingym-
                        dms/TRPC_THEMA_Chan_2017', 'proteingym-
                        dms/UBC9_HUMAN_Weile_2017', 'proteingym-
                        dms/UBE4B_HUMAN_Tsuboyama_2023_3L1X', 'proteingym-
                        dms/UBE4B_HUMAN_Tsuboyama_2023_3L1X_indels',
                        'proteingym-dms/UBE4B_MOUSE_Starita_2013',
                        'proteingym-dms/UBR5_HUMAN_Tsuboyama_2023_1I2T',
                        'proteingym-
                        dms/UBR5_HUMAN_Tsuboyama_2023_1I2T_indels',
                        'proteingym-dms/VG08_BPP22_Tsuboyama_2023_2GP8',
                        'proteingym-
                        dms/VG08_BPP22_Tsuboyama_2023_2GP8_indels',
                        'proteingym-dms/VILI_CHICK_Tsuboyama_2023_1YU5',
                        'proteingym-
                        dms/VILI_CHICK_Tsuboyama_2023_1YU5_indels',
                        'proteingym-dms/VKOR1_HUMAN_Chiasson_2020_abundance',
                        'proteingym-dms/VKOR1_HUMAN_Chiasson_2020_activity',
                        'proteingym-dms/VRPI_BPT7_Tsuboyama_2023_2WNM',
                        'proteingym-dms/VRPI_BPT7_Tsuboyama_2023_2WNM_indels',
                        'proteingym-dms/YAIA_ECOLI_Tsuboyama_2023_2KVT',
                        'proteingym-dms/YAP1_HUMAN_Araya_2012', 'proteingym-
                        dms/YNZC_BACSU_Tsuboyama_2023_2JVD', 'proteingym-
                        dms/YNZC_BACSU_Tsuboyama_2023_2JVD_indels']
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
