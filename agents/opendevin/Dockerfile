FROM biomlbench-env

# where to put submission, will be extracted
ARG SUBMISSION_DIR
ENV SUBMISSION_DIR=${SUBMISSION_DIR}
# where to put any logs, will be extracted
ARG LOGS_DIR
ENV LOGS_DIR=${LOGS_DIR}
# where to put any code, will be extracted
ARG CODE_DIR
ENV CODE_DIR=${CODE_DIR}
# where to put any other agent-specific files, will not be necessarily extracted
ARG AGENT_DIR
ENV AGENT_DIR=${AGENT_DIR}

RUN mkdir -p ${LOGS_DIR} ${AGENT_DIR} ${CODE_DIR} ${SUBMISSION_DIR}

ARG CONDA_ENV_NAME=agent
ENV DEBIAN_FRONTEND=noninteractive

# Install Docker and system dependencies needed for OpenHands
RUN curl -fsSL https://get.docker.com -o /tmp/get-docker.sh && \
    chmod 700 /tmp/get-docker.sh && \
    /tmp/get-docker.sh && \
    sudo usermod -aG docker nonroot

# Install additional system packages that OpenHands might need
RUN apt-get update && apt-get install -y \
    git \
    curl \
    build-essential \
    nodejs \
    npm \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python dependencies via pip (avoid conda TOS issues)
COPY requirements.txt ${AGENT_DIR}/requirements.txt
RUN conda run -n ${CONDA_ENV_NAME} pip install --upgrade pip && \
    conda run -n ${CONDA_ENV_NAME} pip install -r ${AGENT_DIR}/requirements.txt --upgrade && \
    conda clean -afy

# Copy the local OpenHands code (but don't install it as a package to avoid dependency conflicts)
COPY ./ref/OpenHands ${AGENT_DIR}/OpenHands

# Copy OpenDevin agent files (including additional_notes.txt)
COPY setup.py start.py templates.py utils.py build.sh additional_notes.txt ${AGENT_DIR}/
COPY start.sh ${AGENT_DIR}/
COPY entrypoint.sh /agent_entrypoint.sh

RUN chmod +x /agent_entrypoint.sh ${AGENT_DIR}/build.sh ${AGENT_DIR}/start.sh

# Configure passwordless sudo for the build script and conda commands
RUN echo "nonroot ALL=(ALL) NOPASSWD: ${AGENT_DIR}/build.sh" >> /etc/sudoers && \
    echo "ALL ALL=(ALL) NOPASSWD: ${AGENT_DIR}/build.sh" >> /etc/sudoers && \
    echo "nonroot ALL=(ALL) NOPASSWD: /opt/conda/bin/conda" >> /etc/sudoers

# run the agent's entrypoint script instead (which also runs root entrypoint)
ENTRYPOINT ["/agent_entrypoint.sh"]
