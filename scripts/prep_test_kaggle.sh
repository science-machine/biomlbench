#!/bin/bash

# Script to prepare and test all Kaggle tasks in biomlbench
# This script will:
# 1. Prepare each task
# 2. Run the dummy agent on each task  
# 3. Grade the results
# 4. Collect summary statistics

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Create results directory with timestamp
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
RESULTS_DIR="$PROJECT_ROOT/results/kaggle_batch_$TIMESTAMP"
mkdir -p "$RESULTS_DIR"

# Log file
LOG_FILE="$RESULTS_DIR/batch_test.log"

echo -e "${BLUE}=== BioML Bench Kaggle Tasks Batch Testing ===${NC}"
echo "Results will be saved to: $RESULTS_DIR"
echo "Log file: $LOG_FILE"
echo

# Initialize counters
TOTAL_TASKS=0
SUCCESSFUL_TASKS=0
FAILED_TASKS=0
declare -a FAILED_TASK_LIST=()

# Function to log messages
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Function to run a single task
run_task() {
    local task_name="$1"
    local task_id="kaggle/$task_name"
    
    echo -e "${YELLOW}[$((TOTAL_TASKS + 1))] Processing task: $task_id${NC}"
    log "Starting task: $task_id"
    
    # Step 1: Prepare the task
    echo "  → Preparing task..."
    if ! uv run biomlbench prepare -t "$task_id" >> "$LOG_FILE" 2>&1; then
        echo -e "  ${RED}✗ Failed to prepare task${NC}"
        log "ERROR: Failed to prepare task $task_id"
        return 1
    fi
    
    # Step 2: Run dummy agent
    echo "  → Running dummy agent..."
    local run_output
    run_output=$(uv run biomlbench run-agent --agent dummy --task-id "$task_id" 2>&1)
    echo "$run_output" >> "$LOG_FILE"
    
    # Extract submission file path from output
    local submission_file
    submission_file=$(echo "$run_output" | grep "Submission file ready for grading:" | sed 's/.*: //')
    
    if [ -z "$submission_file" ] || [ ! -f "$submission_file" ]; then
        echo -e "  ${RED}✗ Failed to run dummy agent or generate submission${NC}"
        log "ERROR: Failed to run dummy agent for task $task_id"
        return 1
    fi
    
    # Step 3: Grade the results
    echo "  → Grading results..."
    local grade_output_dir="$RESULTS_DIR/$task_name"
    mkdir -p "$grade_output_dir"
    
    if ! uv run biomlbench grade --submission "$submission_file" --output-dir "$grade_output_dir" >> "$LOG_FILE" 2>&1; then
        echo -e "  ${RED}✗ Failed to grade results${NC}"
        log "ERROR: Failed to grade results for task $task_id"
        return 1
    fi
    
    echo -e "  ${GREEN}✓ Task completed successfully${NC}"
    log "SUCCESS: Task $task_id completed successfully"
    return 0
}

# Get all task directories (exclude __pycache__ and files)
KAGGLE_TASKS_DIR="$PROJECT_ROOT/biomlbench/tasks/kaggle"
cd "$KAGGLE_TASKS_DIR"

echo "Scanning for tasks in: $KAGGLE_TASKS_DIR"
echo

# Process each task directory
for task_dir in */; do
    # Remove trailing slash and skip __pycache__
    task_name="${task_dir%/}"
    
    # Skip non-task directories
    if [ "$task_name" = "__pycache__" ] || [ ! -d "$task_name" ]; then
        continue
    fi
    
    TOTAL_TASKS=$((TOTAL_TASKS + 1))
    
    if run_task "$task_name"; then
        SUCCESSFUL_TASKS=$((SUCCESSFUL_TASKS + 1))
    else
        FAILED_TASKS=$((FAILED_TASKS + 1))
        FAILED_TASK_LIST+=("$task_name")
    fi
    
    echo
done

# Generate summary report
echo -e "${BLUE}=== BATCH TEST SUMMARY ===${NC}"
echo "Total tasks processed: $TOTAL_TASKS"
echo -e "Successful: ${GREEN}$SUCCESSFUL_TASKS${NC}"
echo -e "Failed: ${RED}$FAILED_TASKS${NC}"

if [ ${#FAILED_TASK_LIST[@]} -gt 0 ]; then
    echo -e "\n${RED}Failed tasks:${NC}"
    for failed_task in "${FAILED_TASK_LIST[@]}"; do
        echo "  - $failed_task"
    done
fi

# Save summary to file
SUMMARY_FILE="$RESULTS_DIR/batch_summary.txt"
{
    echo "BioML Bench Kaggle Tasks Batch Test Summary"
    echo "Timestamp: $TIMESTAMP"
    echo "Total tasks: $TOTAL_TASKS"
    echo "Successful: $SUCCESSFUL_TASKS" 
    echo "Failed: $FAILED_TASKS"
    echo "Success rate: $(( SUCCESSFUL_TASKS * 100 / TOTAL_TASKS ))%"
    echo
    if [ ${#FAILED_TASK_LIST[@]} -gt 0 ]; then
        echo "Failed tasks:"
        for failed_task in "${FAILED_TASK_LIST[@]}"; do
            echo "  - $failed_task"
        done
    fi
} > "$SUMMARY_FILE"

echo -e "\nResults saved to: ${BLUE}$RESULTS_DIR${NC}"
echo -e "Summary: ${BLUE}$SUMMARY_FILE${NC}"
echo -e "Log file: ${BLUE}$LOG_FILE${NC}"

log "Batch testing completed. Success rate: $(( SUCCESSFUL_TASKS * 100 / TOTAL_TASKS ))%"

# Exit with error code if any tasks failed
if [ $FAILED_TASKS -gt 0 ]; then
    exit 1
fi

echo -e "\n${GREEN}✓ All tasks completed successfully!${NC}" 