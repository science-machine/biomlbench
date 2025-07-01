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
                        expression-prediction-dataset-task']
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
                        expression-prediction-dataset-task']

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
                        expression-prediction-dataset-task']
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
                        expression-prediction-dataset-task']
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
