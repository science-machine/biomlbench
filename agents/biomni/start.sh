#!/bin/bash
set -x # Print commands and their arguments as they are executed

cd ${AGENT_DIR}

eval "$(conda shell.bash hook)" # make conda available to the shell
conda activate biomni_e1  

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
cp /home/instructions.txt ${AGENT_DIR}/instructions.txt

# add agent-specific instructions with a linebreak in between
echo "" >> ${AGENT_DIR}/instructions.txt
envsubst < ${AGENT_DIR}/additional_notes.txt >> ${AGENT_DIR}/instructions.txt

# append the task instructions with a linebreak in between
printf "\nCOMPETITION INSTRUCTIONS\n------\n\n" >> ${AGENT_DIR}/instructions.txt

# overwrite description.md with description_obfuscated.md if $OBFUSCATE is set
if [ "$OBFUSCATE" = "true" ]; then
  if [ ! -w /home/data/ ]; then
    echo "Obfuscation not implemented for read-only mounts"
    exit 1
  fi
  mv /home/data/description_obfuscated.md /home/data/description.md
fi
cat /home/data/description.md >> ${AGENT_DIR}/instructions.txt

# Create workspace and logs directories
mkdir -p ${AGENT_DIR}/logs
mkdir -p ${AGENT_DIR}/workspaces

# Verify that data lake files are present
echo "Verifying pre-downloaded data lake files..."
python ${AGENT_DIR}/verify_data_download.py
if [ $? -ne 0 ]; then
    echo "Warning: Some data lake files may be missing. Biomni will download them at runtime."
fi

# run with timeout, using the Biomni format with correct parameters
timeout $TIME_LIMIT_SECS python main.py \
  --path ${AGENT_DIR} \
  "$@"

if [ $? -eq 124 ]; then
  echo "Timed out after $TIME_LIMIT"
fi

# Copy results to expected output directories
echo "Copying Biomni results..."

# Copy logs
if [ -d "${AGENT_DIR}/logs" ]; then
  cp -r ${AGENT_DIR}/logs/* ${LOGS_DIR}/ 2>/dev/null || true
fi

# Also check the default logs directory that Biomni creates
if [ -d "./logs" ]; then
  cp -r ./logs/* ${LOGS_DIR}/ 2>/dev/null || true
fi

# Copy any generated code files
if [ -d "${AGENT_DIR}/workspaces" ]; then
  find ${AGENT_DIR}/workspaces -name "*.py" -type f -exec cp {} ${CODE_DIR}/ \; 2>/dev/null || true
fi

echo "Biomni execution complete."