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

RUN mkdir ${LOGS_DIR} ${CODE_DIR} ${AGENT_DIR}

ARG CONDA_ENV_NAME=agent

# Accept Anaconda TOS for defaults channels (needed for non-interactive builds)
SHELL ["/bin/bash", "-lc"]
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Clone the Biomni repo and create the environment
# TODO: Once they have a stable release, switch to using that
RUN git clone https://github.com/snap-stanford/Biomni.git ${AGENT_DIR}/Biomni && \
    cd ${AGENT_DIR}/Biomni && \
    git checkout 0f9dbf925a0512292a2486c87f87152484a9d18a && \
    cd ${AGENT_DIR}/Biomni/biomni_env && \
    conda env create -f environment.yml && \
    cd .. 

RUN conda run -n biomni_e1 pip install ${AGENT_DIR}/Biomni

COPY pre_download_data.py ${AGENT_DIR}/pre_download_data.py

RUN echo "🚀 Starting Biomni data lake pre-download..." && \
    echo "�� This will download ~80+ biomedical datasets (~20 GB total)" && \
    echo "⏱️  This step takes around an hour..." && \
    conda run -n biomni_e1 python ${AGENT_DIR}/pre_download_data.py && \
    mkdir -p ${AGENT_DIR}/biomni_data && \
    cp -r /tmp/biomni_data/* ${AGENT_DIR}/biomni_data/ && \
    rm -rf /tmp/biomni_data && \
    echo "✅ Data lake pre-download completed successfully!"

# Fix NumPy version compatibility with RDKit AFTER data lake download
RUN conda run -n biomni_e1 pip install "numpy<2.0" --force-reinstall && \
    echo "✅ NumPy downgraded for RDKit compatibility"

COPY . ${AGENT_DIR}


