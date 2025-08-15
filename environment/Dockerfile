FROM ubuntu:22.04

# Avoid interactive dialog from apt-get and other packages requiring configuration
ENV DEBIAN_FRONTEND=noninteractive

# install basic packages with updated versions
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    git \
    vim \
    nano \
    unzip \
    zip \
    p7zip-full \
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    python-is-python3 \
    build-essential \
    openssh-server \
    tmux \
    gettext \
    sudo \
    ffmpeg \
    libsm6 \
    libxext6 \
    libstdc++6 \
    gcc-9 \
    g++-9 \
    && pip install jupyter \
    && rm -rf /var/lib/apt/lists/* # removes cache

RUN pip install virtualenv \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh \
    && /opt/conda/bin/conda init

ARG CONDA_ENV_NAME=agent
ARG PYTHON_VERSION=3.11
ARG REQUIREMENTS=/tmp/requirements.txt

COPY environment/requirements.txt ${REQUIREMENTS}

# create conda environment and optionally install the requirements to it
ENV CONDA_PLUGINS_AUTO_ACCEPT_TOS=true
RUN /opt/conda/bin/conda create -n ${CONDA_ENV_NAME} python=${PYTHON_VERSION} -y
ARG INSTALL_HEAVY_DEPENDENCIES=true
ENV INSTALL_HEAVY_DEPENDENCIES=${INSTALL_HEAVY_DEPENDENCIES}

# Install biomedical dependencies with correct GCC version for RDKit
RUN if [ "$INSTALL_HEAVY_DEPENDENCIES" = "true" ]; then \
    /opt/conda/bin/conda install -n ${CONDA_ENV_NAME} -c conda-forge gcc=12.1.0 -y && \
    /opt/conda/bin/conda install -n ${CONDA_ENV_NAME} -c conda-forge rdkit biopython scikit-learn pandas numpy matplotlib seaborn -y && \
    /opt/conda/bin/conda run -n ${CONDA_ENV_NAME} pip install -r /tmp/requirements.txt && \
    /opt/conda/bin/conda run -n ${CONDA_ENV_NAME} pip install tensorflow[and-cuda]==2.17 && \
    /opt/conda/bin/conda run -n ${CONDA_ENV_NAME} pip install torch==2.2.0 torchaudio==2.2.0 torchtext==0.17.0 torchvision==0.17.0 && \
    /opt/conda/bin/conda clean -afy ; fi

ENV PATH="/opt/conda/bin:${PATH}"

# Installs from here onward go into the conda base env; previous was installed to /usr/bin/python

# Install stuff for the grading server: biomlbench and flask
# Copy only necessary files for package installation (improves Docker layer caching)
COPY pyproject.toml README.md /biomlbench/
COPY biomlbench /biomlbench/biomlbench
COPY scripts /biomlbench/scripts
RUN /opt/conda/bin/conda create -n biomlb python=3.11 -y
RUN /opt/conda/bin/conda run -n biomlb pip install flask && \
    /opt/conda/bin/conda run -n biomlb pip install -e /biomlbench

# Reset DEBIAN_FRONTEND
ENV DEBIAN_FRONTEND=

# Make private directory (root) owner-only. Grading server will be added here, later in the build
# The test set answers will be added here separately via a mounted docker volume
RUN mkdir /private && chmod 700 /private

# Copy over relevant files
COPY environment/grading_server.py /private/grading_server.py
COPY environment/instructions.txt /home/instructions.txt
COPY environment/instructions_obfuscated.txt /home/instructions_obfuscated.txt
COPY environment/validate_submission.sh /home/validate_submission.sh
COPY environment/entrypoint.sh /entrypoint.sh

# Create nonroot user; make entrypoint executable
RUN useradd -m nonroot \
    && mkdir /home/submission \
    && chown nonroot:nonroot /home/submission \
    && chmod 755 /home/submission \
    && chmod +x /entrypoint.sh

WORKDIR /home

# IMPORTANT: This needs to run as root! Downstream Dockerfiles must not change the default USER for when the container starts.
# Entrypoint script is in charge of setting up the user environment and running the grading server
ENTRYPOINT ["/entrypoint.sh"]
