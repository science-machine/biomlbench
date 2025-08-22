#!/bin/bash

# BioML-bench Task Preparation Wrapper Script
# 
# This is a convenient wrapper around prepare_all_tasks.py that provides
# some common usage patterns and environment setup.

set -e  # Exit on any error

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Change to project root directory
cd "$PROJECT_ROOT"

echo "BioML-bench Task Preparation Script"
echo "==================================="
echo "Project root: $PROJECT_ROOT"
echo

# Default values
THREADS=4
JOBS_FILE="scripts/gcp-deploy/production-jobs.txt"

# Parse command line arguments for common patterns
case "${1:-}" in
    "help"|"--help"|"-h")
        echo "Usage: $0 [COMMAND] [OPTIONS]"
        echo
        echo "Commands:"
        echo "  all [threads]     - Prepare all tasks (default: 4 threads)"
        echo "  quick [threads]   - Quick preparation (skip verification, 8 threads)"
        echo "  retry [threads]   - Retry failed tasks from failed_tasks.txt"
        echo "  test              - Test with a small subset (3 tasks, 2 threads)"
        echo "  help              - Show this help"
        echo
        echo "Examples:"
        echo "  $0 all 8          - Prepare all tasks with 8 threads"
        echo "  $0 quick          - Quick preparation with 8 threads"
        echo "  $0 retry 2        - Retry failed tasks with 2 threads"
        echo "  $0 test           - Test run with 3 tasks"
        echo
        echo "For advanced usage, run the Python script directly:"
        echo "  python scripts/prepare_all_tasks.py --help"
        exit 0
        ;;
    "all")
        THREADS=${2:-4}
        echo "Preparing ALL tasks with $THREADS threads..."
        python scripts/prepare_all_tasks.py --threads "$THREADS"
        ;;
    "quick")
        THREADS=${2:-8}
        echo "Quick preparation (skip verification) with $THREADS threads..."
        python scripts/prepare_all_tasks.py --threads "$THREADS" --skip-verification
        ;;
    "retry")
        THREADS=${2:-2}
        if [[ ! -f "failed_tasks.txt" ]]; then
            echo "Error: failed_tasks.txt not found. Run a preparation first."
            exit 1
        fi
        echo "Retrying failed tasks with $THREADS threads..."
        python scripts/prepare_all_tasks.py --task-list failed_tasks.txt --threads "$THREADS"
        ;;
    "test")
        echo "Test run: preparing 3 tasks with 2 threads..."
        # Create a small test file with just a few tasks
        head -20 "$JOBS_FILE" | grep -v '^#' | head -6 > test_tasks_temp.txt
        python scripts/prepare_all_tasks.py --jobs test_tasks_temp.txt --threads 2
        rm -f test_tasks_temp.txt
        ;;
    "")
        # Default behavior - prepare all with default threads
        echo "Preparing ALL tasks with $THREADS threads (default)..."
        echo "Use '$0 help' for more options."
        python scripts/prepare_all_tasks.py --threads "$THREADS"
        ;;
    *)
        echo "Unknown command: $1"
        echo "Use '$0 help' for available commands."
        exit 1
        ;;
esac

echo
echo "Task preparation completed!"
echo "Check task_preparation.log for detailed logs."

# Show summary if failed_tasks.txt was created
if [[ -f "failed_tasks.txt" ]]; then
    FAILED_COUNT=$(wc -l < failed_tasks.txt)
    if [[ $FAILED_COUNT -gt 0 ]]; then
        echo
        echo "Warning: $FAILED_COUNT tasks failed. Retry with:"
        echo "  $0 retry"
    fi
fi 