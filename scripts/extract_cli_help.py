#!/usr/bin/env python3
"""
Extract real CLI help output for documentation.

This script runs the actual CLI commands with --help to capture the real
argparse help output and formats it as markdown documentation.
"""

import subprocess
import sys
from pathlib import Path


def run_help_command(cmd):
    """Run a command with --help and capture output."""
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        if result.returncode == 0 or result.stderr == "":
            return result.stdout
        else:
            return result.stderr if result.stderr else result.stdout
    except subprocess.TimeoutExpired:
        return "Command timed out"
    except Exception as e:
        return f"Error running command: {e}"


def generate_comprehensive_cli_docs():
    """Generate comprehensive CLI documentation with real help output."""

    md_content = """# CLI Reference

BioML-bench's command-line interface for managing biomedical benchmark tasks and running agent evaluations.

::: biomlbench.cli
    options:
      show_source: false
      heading_level: 2

## Main Help

```
"""

    # Get main help
    main_help = run_help_command([sys.executable, "-m", "biomlbench.cli", "--help"])
    md_content += main_help + "\n```\n\n"

    # Get help for each subcommand
    commands = ["prepare", "run-agent", "grade", "grade-sample", "run-baseline", "dev"]

    md_content += "## Commands\n\n"

    for cmd in commands:
        md_content += f"### `biomlbench {cmd}`\n\n"

        help_output = run_help_command([sys.executable, "-m", "biomlbench.cli", cmd, "--help"])
        md_content += f"```\n{help_output}\n```\n\n"

        # Special handling for dev subcommands
        if cmd == "dev":
            dev_commands = ["download-leaderboard", "prepare-human-baselines"]
            for dev_cmd in dev_commands:
                md_content += f"#### `biomlbench dev {dev_cmd}`\n\n"
                dev_help = run_help_command(
                    [sys.executable, "-m", "biomlbench.cli", "dev", dev_cmd, "--help"]
                )
                md_content += f"```\n{dev_help}\n```\n\n"

    # Add practical examples
    md_content += """## Usage Examples

### Task Preparation

```bash
# Prepare a specific task
biomlbench prepare -t caco2-wang

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
biomlbench run-agent --agent dummy --task-id caco2-wang

# Run agent on multiple tasks with parallel workers
biomlbench run-agent \\
    --agent aide \\
    --task-list experiments/splits/caco2-wang.txt \\
    --n-workers 4 \\
    --n-seeds 3

# Run with custom container configuration
biomlbench run-agent \\
    --agent aide \\
    --task-id caco2-wang \\
    --container-config custom_config.json \\
    --retain-container
```

### Evaluation and Grading

```bash
# Grade multiple task submissions
biomlbench grade \\
    --submission runs/my-run-group/submission.jsonl \\
    --output-dir results/

# Grade single task submission
biomlbench grade-sample submission.csv caco2-wang

# Run and grade baselines
biomlbench run-baseline caco2-wang --baseline all
biomlbench grade \\
    --submission baseline_submissions/submission.jsonl \\
    --output-dir baseline_results/
```

### Development Commands

```bash
# Download leaderboards for all tasks
biomlbench dev download-leaderboard --all

# Prepare human baselines for a specific task
biomlbench dev prepare-human-baselines -t caco2-wang --force
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
biomlbench run-agent --agent dummy --task-list experiments/splits/caco2-wang.txt

# 3. Grade results (submission.jsonl is auto-generated)
biomlbench grade \\
    --submission runs/latest-run-group/submission.jsonl \\
    --output-dir results/
```

### Debugging Agent Issues

```bash
# Run with container retention for debugging
biomlbench run-agent \\
    --agent my-agent \\
    --task-id caco2-wang \\
    --retain-container

# Check logs in the run directory
ls runs/latest-run-group/*/logs/
```
"""

    return md_content


def main():
    """Generate comprehensive CLI documentation."""

    # Generate the documentation
    md_content = generate_comprehensive_cli_docs()

    # Write to file
    docs_file = Path(__file__).parent.parent / "docs" / "api" / "cli.md"
    with open(docs_file, "w") as f:
        f.write(md_content)

    print(f"Generated CLI documentation: {docs_file}")


if __name__ == "__main__":
    main()
