# Failed Runs Analysis Tools

This directory contains tools for analyzing failed biomlbench runs.

## download_failed_runs.py

Downloads and extracts all failed run artifacts from the S3 bucket for local analysis.

### Usage

```bash
# Basic usage - download to ./failed_runs directory
python scripts/analyze_failures/download_failed_runs.py

# Download to a custom directory
python scripts/analyze_failures/download_failed_runs.py --output-dir /path/to/analysis

# Dry run to see what would be downloaded
python scripts/analyze_failures/download_failed_runs.py --dry-run

# Keep compressed files after extraction
python scripts/analyze_failures/download_failed_runs.py --keep-compressed

# Use more concurrent downloads (default: 5)
python scripts/analyze_failures/download_failed_runs.py --concurrent 10
```

### Features

- **Parallel Downloads**: Downloads multiple files concurrently for speed
- **Automatic Extraction**: Extracts tar.gz files automatically
- **Progress Tracking**: Shows detailed progress and success/failure rates  
- **Organized Output**: Creates separate directories for each failed run
- **Cleanup**: Removes compressed files after extraction (optional)
- **Dry Run**: Preview what would be downloaded without actually doing it

### Output Structure

The script creates a directory structure like:
```
failed_runs/
├── run-group_agent1_12345678/
│   ├── logs/
│   ├── outputs/
│   └── ...
├── run-group_agent2_87654321/
│   ├── logs/
│   ├── outputs/
│   └── ...
└── ...
```

Each directory contains the complete failed run artifacts for analysis.

### Requirements

- AWS CLI configured with access to the biomlbench S3 bucket
- `tar` command available (standard on Linux/macOS)
- Python 3.6+ with standard library modules

## parse_aide_journals.py

Parses AIDE agent journal files from failed runs and formats them into readable plain text files.

### Usage

```bash
# Basic usage - parse all AIDE journals from failed_runs directory
python scripts/analyze_failures/parse_aide_journals.py

# Parse from custom directory
python scripts/analyze_failures/parse_aide_journals.py --failed-runs-dir /path/to/failed_runs

# Dry run to see what would be processed
python scripts/analyze_failures/parse_aide_journals.py --dry-run

# Use more concurrent processing
python scripts/analyze_failures/parse_aide_journals.py --concurrent 10
```

### Features

- **Automatic Discovery**: Finds all AIDE journal.json files in the failed runs
- **Plain Text Output**: Creates clean, readable text files for human review
- **Organized Structure**: Groups by run ID and creates separate files per task
- **Parallel Processing**: Processes multiple journals concurrently
- **Comprehensive Parsing**: Extracts plans, code, outputs, errors, and analysis from each step

### Output Structure

The script creates organized output like:
```
aide_journals/
├── 2025-08-21T02-06-57-GMT_run-group_aide/
│   └── kaggle_uw-madison-gi-tract-image-segmentation.txt
├── 2025-08-20T17-09-05-GMT_run-group_aide/
│   └── kaggle_uw-madison-gi-tract-image-segmentation.txt
├── 2025-08-20T16-37-58-GMT_run-group_aide/
│   └── manual_open-problems-label-projection.txt
└── ...
```

### What It Extracts

For each step in the agent's execution, the parser extracts:
- **Plan**: The agent's reasoning and plan for the step
- **Code**: The actual code being executed  
- **Output**: Terminal output and results
- **Analysis**: Agent's analysis of the results
- **Errors**: Exception types and details when failures occur
- **Timing**: Execution timestamps and durations
- **Status**: Whether the step was marked as buggy

### Text File Format

Each text file includes:
- Header with run metadata (run ID, task, timestamps, etc.)
- Each step clearly separated with dividers
- Clean formatting with emojis for easy scanning
- All relevant data for understanding what went wrong

### Requirements

- Python 3.6+ with standard library modules
- Access to the failed_runs directory with AIDE journal files 