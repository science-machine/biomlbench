#!/bin/bash
set -x # Print commands and their arguments as they are executed

# Navigate to the agent directory
cd ${AGENT_DIR}

# Activate the conda environment
eval "$(conda shell.bash hook)"
conda activate autobioml

# Set PYTHONPATH to include the AutoBioML source
export PYTHONPATH="${AGENT_DIR}/AgBioML-Challenge/src:${PYTHONPATH}"

# Parse command line arguments
MODEL="gpt-4o"
MAX_ITERATIONS=25
ENABLE_PUBLIC_EVAL="true"

while [[ $# -gt 0 ]]; do
    case $1 in
        --model)
            MODEL="$2"
            shift 2
            ;;
        --max_iterations)
            MAX_ITERATIONS="$2"
            shift 2
            ;;
        --enable_public_evaluation)
            ENABLE_PUBLIC_EVAL="$2"
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done

# Create a working directory for the task
WORK_DIR="/home/autobioml_work"
mkdir -p ${WORK_DIR}

# Run the adapter to set up the task
python ${AGENT_DIR}/adapter.py \
    --task-id "${TASK_ID}" \
    --data-dir "/home/data" \
    --work-dir "${WORK_DIR}" \
    --submission-dir "${SUBMISSION_DIR}"

# Check if adapter succeeded
if [ $? -ne 0 ]; then
    echo "Error: Failed to adapt task for AutoBioML"
    exit 1
fi

# Change to the working directory
cd ${WORK_DIR}

# Run AutoBioML using our wrapper that handles local execution
# The adapter creates a challenges/current directory with the task
python ${AGENT_DIR}/autobioml_wrapper.py \
    task_dir=challenges/current \
    model=${MODEL} \
    max_iterations=${MAX_ITERATIONS} \
    enable_public_evaluation=${ENABLE_PUBLIC_EVAL} \
    docker.enabled=false \
    hydra.run.dir="${WORK_DIR}/outputs"

# Check if AutoBioML succeeded
if [ $? -ne 0 ]; then
    echo "Error: AutoBioML execution failed"
    exit 1
fi

# Find the latest output directory
LATEST_OUTPUT=$(find ${WORK_DIR}/outputs -type d -name "[0-9]*" | sort -r | head -n 1)

if [ -z "$LATEST_OUTPUT" ]; then
    echo "Error: No output directory found"
    exit 1
fi

echo "AutoBioML output directory: ${LATEST_OUTPUT}"

# Convert AutoBioML output to BioML-bench submission format
python ${AGENT_DIR}/convert_output.py \
    --autobioml-output "${LATEST_OUTPUT}" \
    --submission-dir "${SUBMISSION_DIR}" \
    --task-id "${TASK_ID}"

# Copy logs
if [ -f "${LATEST_OUTPUT}/lab_notebook.md" ]; then
    cp "${LATEST_OUTPUT}/lab_notebook.md" "${LOGS_DIR}/"
fi

# Copy any generated code
if [ -d "${LATEST_OUTPUT}" ]; then
    find "${LATEST_OUTPUT}" -name "*.py" -o -name "*.ipynb" | while read file; do
        cp "$file" "${CODE_DIR}/" 2>/dev/null || true
    done
fi

echo "AutoBioML agent execution completed" 