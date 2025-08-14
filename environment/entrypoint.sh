#!/bin/bash

# Print commands and their arguments as they are executed
set -x

{
  # log into /home/logs
  LOGS_DIR=/home/logs
  mkdir -p $LOGS_DIR

  echo "Starting entrypoint.sh at $(date)"
  echo "Directory contents before chmod:"
  ls -la /home/

  # chmod the /home directory such that nonroot users can work on everything within it. We do this at container start
  # time so that anything added later in agent-specific Dockerfiles will also receive the correct permissions.
  # (this command does `chmod a+rw /home` but with the exception of /home/data, which is a read-only volume)
  # Add timeout to prevent hanging on large directory trees
  echo "Starting chmod operation at $(date)"
  timeout 60s find /home -path /home/data -prune -o -exec chmod a+rw {} \; || {
    echo "chmod operation timed out or failed at $(date)"
    echo "This is likely due to a large conda environment. Continuing anyway..."
  }
  echo "chmod operation completed at $(date)"
  
  ls -l /home

  echo "Starting grading server at $(date)"
  # Launch grading server, stays alive throughout container lifetime to service agent requests.
  /opt/conda/bin/conda run -n biomlb python /private/grading_server.py
} 2>&1 | tee $LOGS_DIR/entrypoint.log
