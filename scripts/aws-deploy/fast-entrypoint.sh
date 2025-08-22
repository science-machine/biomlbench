#!/bin/bash

# Print commands and their arguments as they are executed
set -x

{
  # log into /home/logs
  LOGS_DIR=/home/logs
  mkdir -p $LOGS_DIR

  echo "Starting FAST entrypoint.sh at $(date)"
  echo "Directory contents before setup:"
  ls -la /home/

  # FAST permission setup - only what's absolutely necessary
  echo "Starting minimal permission setup at $(date)"
  
  # Make /home accessible
  chmod 755 /home 2>/dev/null || true
  
  # Only set up directories that actually exist and are needed
  if [ -d "/home/runner" ]; then
    chown -R nonroot:nonroot /home/runner 2>/dev/null || true
  fi
  
  if [ -d "/home/agent" ]; then
    chown -R nonroot:nonroot /home/agent 2>/dev/null || true
    chmod -R u+rw /home/agent 2>/dev/null || true
  fi
  
  if [ -d "/home/logs" ]; then
    chown nonroot:nonroot /home/logs 2>/dev/null || true
  fi
  
  if [ -d "/home/submission" ]; then
    chown nonroot:nonroot /home/submission 2>/dev/null || true
  fi
  
  if [ -d "/home/code" ]; then
    chown -R nonroot:nonroot /home/code 2>/dev/null || true
    chmod -R u+rw /home/code 2>/dev/null || true
  fi
  
  # Make sure nonroot can work in /home
  chmod g+w /home 2>/dev/null || true
  
  echo "Minimal permission setup completed at $(date)"
  
  ls -l /home

  echo "Starting grading server at $(date)"
  # Launch grading server, stays alive throughout container lifetime to service agent requests.
  /opt/conda/bin/conda run -n biomlb python /private/grading_server.py
} 2>&1 | tee $LOGS_DIR/entrypoint.log 