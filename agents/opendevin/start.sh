#!/bin/bash

# Print commands and their arguments as they are executed
set -x

cd ${AGENT_DIR}

eval "$(conda shell.bash hook)" # make conda available to the shell
conda activate agent

# determine hardware available
if command -v nvidia-smi &> /dev/null && nvidia-smi --query-gpu=name --format=csv,noheader &> /dev/null; then
  HARDWARE=$(nvidia-smi --query-gpu=name --format=csv,noheader \
    | sed 's/^[ \t]*//' \
    | sed 's/[ \t]*$//' \
    | sort \
    | uniq -c \
    | sed 's/^ *\([0-9]*\) *\(.*\)$/\1 \2/' \
    | paste -sd ', ' -)
else
  HARDWARE="a CPU"
fi
export HARDWARE

# check that we can use the GPU in PyTorch
python -c "import torch; print(torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'WARNING: No GPU')"
# check that we can use the GPU in TensorFlow
python -c "import tensorflow as tf; print('GPUs Available: ', tf.config.list_physical_devices('GPU'))"

# convert $TIME_LIMIT_SECS to more readable format for prompt
format_time() {
  local time_in_sec=$1
  local hours=$((time_in_sec / 3600))
  local minutes=$(((time_in_sec % 3600) / 60))
  local seconds=$((time_in_sec % 60))
  echo "${hours}hrs ${minutes}mins ${seconds}secs"
}
export TIME_LIMIT=$(format_time $TIME_LIMIT_SECS)

# overwrite instructions.txt with instructions_obfuscated.txt if $OBFUSCATE is set
if [ "$OBFUSCATE" = "true" ]; then
  if [ ! -w /home/data/ ]; then
    echo "Obfuscation not implemented for read-only mounts"
    exit 1
  fi
  mv /home/instructions_obfuscated.txt /home/instructions.txt
fi

# start a new file to store the full instructions, starting with general instructions
cp /home/instructions.txt ${AGENT_DIR}/full_instructions.txt

# add agent-specific instructions with a linebreak in between
echo "" >> ${AGENT_DIR}/full_instructions.txt
envsubst < ${AGENT_DIR}/additional_notes.txt >> ${AGENT_DIR}/full_instructions.txt

# append the task instructions with a linebreak in between
printf "\nCOMPETITION INSTRUCTIONS\n------\n\n" >> ${AGENT_DIR}/full_instructions.txt

# overwrite description.md with description_obfuscated.md if $OBFUSCATE is set
if [ "$OBFUSCATE" = "true" ]; then
  if [ ! -w /home/data/ ]; then
    echo "Obfuscation not implemented for read-only mounts"
    exit 1
  fi
  mv /home/data/description_obfuscated.md /home/data/description.md
fi
cat /home/data/description.md >> ${AGENT_DIR}/full_instructions.txt

# Create workspace and logs directories
mkdir -p ${AGENT_DIR}/logs
mkdir -p ${AGENT_DIR}/workspace

# Set up OpenHands environment
export OPENAI_API_KEY=${OPENAI_API_KEY}
export ANTHROPIC_API_KEY=${ANTHROPIC_API_KEY}

# Create OpenHands config
python3 -c "
import os
import toml
from pathlib import Path

# Get arguments from environment/command line
import sys
args = {}
for arg in sys.argv[1:]:
    if '=' in arg:
        key, value = arg.split('=', 1)
        # Convert dotted keys to nested dict
        keys = key.split('.')
        current = args
        for k in keys[:-1]:
            if k not in current:
                current[k] = {}
            current = current[k]
        current[keys[-1]] = value

# Set default values
config = {
    'core': {
        'workspace_base': '/home/data',
        'workspace_mount_path': '/home/data',
        'workspace_mount_path_in_sandbox': '/home/data',
        'default_agent': args.get('agent', 'CodeActAgent'),
        'max_iterations': int(args.get('max_steps', os.getenv('STEP_LIMIT', 500))),
        'persist_sandbox': False
    },
    'llm': {
        'model': args.get('model', 'gpt-4o'),
        'api_key': os.getenv('OPENAI_API_KEY') or os.getenv('ANTHROPIC_API_KEY'),
        'embedding_model': 'local'
    },
    'sandbox': {
        'user_id': 1000,
        'use_host_network': True,
        'timeout': int(args.get('exec', {}).get('timeout', 3600))
    }
}

# Write config
with open('${AGENT_DIR}/config.toml', 'w') as f:
    toml.dump(config, f)
" "$@"

# Build (no sudo)
./build.sh \
  && conda run -n agent --no-capture-output python setup.py "$@" \
  && conda run -n agent --no-capture-output python start.py "$@"

if [ $? -eq 124 ]; then
  echo "Timed out after $TIME_LIMIT"
fi

# Copy results to expected output directories
echo "Copying OpenHands results..."

# Copy logs
if [ -d "${AGENT_DIR}/logs" ]; then
  cp -r ${AGENT_DIR}/logs/* ${LOGS_DIR}/ 2>/dev/null || true
fi

# Look for any Python files in the workspace that might be solutions
if [ -d "/home/data" ]; then
  find /home/data -name "*.py" -type f | head -5 | while read file; do
    cp "$file" ${CODE_DIR}/ 2>/dev/null || true
  done
  
  # Look for submission files
  find /home/data -name "submission.csv" -type f -exec cp {} ${CODE_DIR}/ \; -quit 2>/dev/null || true
  find /home/data -name "*submission*" -type f | head -3 | while read file; do
    cp "$file" ${CODE_DIR}/ 2>/dev/null || true
  done
fi

echo "OpenDevin execution complete."
