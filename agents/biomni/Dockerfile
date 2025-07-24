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

# Clone the Biomni repo and create the environment
RUN git clone https://github.com/snap-stanford/Biomni.git ${AGENT_DIR}/Biomni && \
    cd ${AGENT_DIR}/Biomni/biomni_env && \
    conda env create -f environment.yml && \
    cd ..

RUN conda activate biomni_e1 && pip install git+https://github.com/snap-stanford/Biomni.git@main
