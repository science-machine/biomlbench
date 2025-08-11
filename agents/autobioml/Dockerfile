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

# Switch back to nonroot user
USER nonroot

# Create conda environment for AutoBioML
ARG CONDA_ENV_NAME=autobioml
RUN conda create -n ${CONDA_ENV_NAME} python=3.10 -y

# Install UV for the nonroot user
USER root
RUN curl -LsSf https://astral.sh/uv/install.sh | sh && \
    mv /root/.cargo/bin/uv /usr/local/bin/ && \
    chmod +x /usr/local/bin/uv
USER nonroot

# Clone the AutoBioML source code from GitHub
RUN git clone https://github.com/AIxBio/AgBioML-Challenge.git ${AGENT_DIR}/AgBioML-Challenge

# Set up the AutoBioML environment
WORKDIR ${AGENT_DIR}/AgBioML-Challenge

# Copy requirements file first
COPY requirements.txt ${AGENT_DIR}/

# Install dependencies using UV
RUN conda run -n ${CONDA_ENV_NAME} uv pip install -e .

# Install additional dependencies from requirements.txt
RUN conda run -n ${CONDA_ENV_NAME} pip install -r ${AGENT_DIR}/requirements.txt

# Create a wrapper script to ensure autobioml is in PATH
RUN echo '#!/bin/bash\nconda run -n autobioml python -m autobioml "$@"' > /home/nonroot/.local/bin/autobioml && \
    chmod +x /home/nonroot/.local/bin/autobioml

ENV PATH="/home/nonroot/.local/bin:${PATH}"

# Copy adapter and wrapper scripts
COPY adapter.py ${AGENT_DIR}/
COPY convert_output.py ${AGENT_DIR}/
COPY autobioml_wrapper.py ${AGENT_DIR}/

# Copy the start script and make scripts executable
COPY start.sh ${AGENT_DIR}/
RUN chmod +x ${AGENT_DIR}/start.sh ${AGENT_DIR}/autobioml_wrapper.py ${AGENT_DIR}/adapter.py ${AGENT_DIR}/convert_output.py

# Set working directory
WORKDIR ${AGENT_DIR} 