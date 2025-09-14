#!/bin/bash
set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to show usage
show_usage() {
    echo "Usage: $0 <task-file>"
    echo ""
    echo "Prepare all tasks listed in a text file for BioML-bench."
    echo ""
    echo "Arguments:"
    echo "  task-file    Path to text file containing task IDs (one per line)"
    echo ""
    echo "Example:"
    echo "  $0 experiments/biomlbench_v0.1a.txt"
    echo ""
    echo "The task file should contain task IDs like:"
    echo "  polarishub/tdcommons-caco2-wang"
    echo "  kaggle/histopathologic-cancer-detection"
    echo "  manual/open-problems-predict-modality"
    echo ""
}

# Check if task file argument was provided
if [[ $# -ne 1 ]]; then
    echo -e "${RED}❌ Error: Task file argument is required${NC}"
    show_usage
    exit 1
fi

TASK_FILE="$1"

# Check if task file exists
if [[ ! -f "$TASK_FILE" ]]; then
    echo -e "${RED}❌ Error: Task file '$TASK_FILE' not found${NC}"
    exit 1
fi

echo -e "${BLUE}📋 BioML-bench Task Preparation${NC}"
echo "==============================="
echo -e "Task file: ${YELLOW}$TASK_FILE${NC}"
echo ""

# Read tasks from file, skip empty lines and comments
mapfile -t TASK_IDS < <(grep -v '^\s*$' "$TASK_FILE" | grep -v '^\s*#')

if [[ ${#TASK_IDS[@]} -eq 0 ]]; then
    echo -e "${RED}❌ Error: No task IDs found in '$TASK_FILE'${NC}"
    echo "Make sure the file contains task IDs (one per line) and is not empty."
    exit 1
fi

echo -e "${YELLOW}Found ${#TASK_IDS[@]} tasks to prepare:${NC}"
for task_id in "${TASK_IDS[@]}"; do
    echo "  - $task_id"
done
echo ""

# Ask for confirmation
read -p "Proceed with preparing all tasks? This may take a long time and download several GB of data. [y/N] " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo -e "${YELLOW}⚠️  Task preparation cancelled by user${NC}"
    exit 0
fi

echo ""
echo -e "${YELLOW}🚀 Starting task preparation...${NC}"
echo ""

# Track success/failure
declare -a PREPARED_TASKS=()
declare -a FAILED_TASKS=()
declare -a SKIPPED_TASKS=()

# Prepare each task
for i in "${!TASK_IDS[@]}"; do
    task_id="${TASK_IDS[$i]}"
    task_num=$((i + 1))
    
    echo -e "${BLUE}[$task_num/${#TASK_IDS[@]}] Preparing: $task_id${NC}"
    echo "----------------------------------------"
    
    # Run biomlbench prepare command
    if biomlbench prepare -t "$task_id"; then
        echo -e "${GREEN}✅ Successfully prepared: $task_id${NC}"
        PREPARED_TASKS+=("$task_id")
    else
        exit_code=$?
        if [[ $exit_code -eq 2 ]]; then
            echo -e "${YELLOW}⏭️  Already prepared: $task_id${NC}"
            SKIPPED_TASKS+=("$task_id")
        else
            echo -e "${RED}❌ Failed to prepare: $task_id${NC}"
            FAILED_TASKS+=("$task_id")
        fi
    fi
    echo ""
done

# Summary
echo -e "${YELLOW}📊 Preparation Summary${NC}"
echo "======================"

if [[ ${#PREPARED_TASKS[@]} -gt 0 ]]; then
    echo -e "${GREEN}Successfully prepared (${#PREPARED_TASKS[@]}):"
    for task in "${PREPARED_TASKS[@]}"; do
        echo "  ✅ $task"
    done
    echo ""
fi

if [[ ${#SKIPPED_TASKS[@]} -gt 0 ]]; then
    echo -e "${YELLOW}Already prepared/skipped (${#SKIPPED_TASKS[@]}):"
    for task in "${SKIPPED_TASKS[@]}"; do
        echo "  ⏭️  $task"
    done
    echo ""
fi

if [[ ${#FAILED_TASKS[@]} -gt 0 ]]; then
    echo -e "${RED}Failed to prepare (${#FAILED_TASKS[@]}):"
    for task in "${FAILED_TASKS[@]}"; do
        echo "  ❌ $task"
    done
    echo ""
    echo -e "${RED}⚠️  Some tasks failed to prepare. Check the logs above for details.${NC}"
    exit 1
else
    echo -e "${GREEN}🎉 All tasks processed successfully!${NC}"
fi

echo ""
echo "Next steps:"
echo "1. Run agents on prepared tasks: biomlbench run-agent --agent <agent> --task-list $TASK_FILE"
echo "2. Or run individual tasks: biomlbench run-agent --agent <agent> --task-id <task-id>"

