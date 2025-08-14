FROM biomlbench-env

# Set up environment variables
ARG SUBMISSION_DIR
ENV SUBMISSION_DIR=${SUBMISSION_DIR}
ARG LOGS_DIR
ENV LOGS_DIR=${LOGS_DIR}
ARG CODE_DIR
ENV CODE_DIR=${CODE_DIR}
ARG AGENT_DIR
ENV AGENT_DIR=${AGENT_DIR}

# Create necessary directories
RUN mkdir -p ${LOGS_DIR} ${CODE_DIR} ${AGENT_DIR}

# Install system dependencies
USER root
RUN apt-get update && apt-get install -y \
    git \
    curl \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install UV package manager
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.cargo/bin:${PATH}"

# Stay as root to create conda environment and set up everything
# This avoids permission issues during entrypoint chmod

# Set conda TOS auto-acceptance as environment variable
ENV CONDA_PLUGINS_AUTO_ACCEPT_TOS=true



# Install UV for system-wide use
RUN curl -LsSf https://astral.sh/uv/install.sh | sh && \
    mv /root/.local/bin/uv /usr/local/bin/ && \
    chmod +x /usr/local/bin/uv

# Copy requirements file first
COPY requirements.txt ${AGENT_DIR}/

# Create conda environment for AutoBioML with core packages that AutoBioML expects
ARG CONDA_ENV_NAME=base
# Create environment with comprehensive package set similar to challenges/Dockerfile
RUN conda install -n ${CONDA_ENV_NAME} python=3.10 && \
    conda install -n ${CONDA_ENV_NAME} -c bioconda -c conda-forge -c pytorch -y \
        scikit-learn pandas seaborn numpy matplotlib scipy 

# Set up the AutoBioML environment
WORKDIR ${AGENT_DIR}/AgBioML-Challenge


# Copy the AutoBioML source code (local development copy)
COPY AgBioML-Challenge ${AGENT_DIR}/AgBioML-Challenge

# Install dependencies using UV
RUN conda run -n ${CONDA_ENV_NAME} uv pip install -e .

# Install additional dependencies from requirements.txt
RUN conda run -n ${CONDA_ENV_NAME} pip install -r ${AGENT_DIR}/requirements.txt

# Create a wrapper script to ensure autobioml is in PATH
# Put the wrapper in /usr/local/bin so it's system-wide and not affected by chmod
RUN echo '#!/bin/bash\nconda run -n autobioml python -m autobioml "$@"' > /usr/local/bin/autobioml && \
    chmod +x /usr/local/bin/autobioml

# Create minimal user directories that won't cause chmod issues
RUN mkdir -p /home/nonroot/.local && \
    chown -R nonroot:nonroot /home/nonroot/.local && \
    chmod -R 755 /home/nonroot/.local

ENV PATH="/usr/local/bin:${PATH}"

# Switch to root to copy and set permissions
USER root

# Copy adapter and wrapper scripts
COPY adapter.py ${AGENT_DIR}/
COPY convert_output.py ${AGENT_DIR}/
COPY autobioml_wrapper.py ${AGENT_DIR}/

# Copy the start script and make scripts executable
COPY start.sh ${AGENT_DIR}/
RUN chmod +x ${AGENT_DIR}/start.sh ${AGENT_DIR}/autobioml_wrapper.py ${AGENT_DIR}/adapter.py ${AGENT_DIR}/convert_output.py && \
    chown -R nonroot:nonroot ${AGENT_DIR}

# Set working directory
WORKDIR ${AGENT_DIR}

# IMPORTANT: Keep running as root for entrypoint.sh 
# The entrypoint will handle user switching and permissions 