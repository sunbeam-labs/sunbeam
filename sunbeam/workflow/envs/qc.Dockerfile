FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/qc_env

COPY sunbeam/workflow/envs/qc.yml ./

# Install environment
RUN mamba env create --file qc.yml --name qc

ENV PATH="/opt/conda/envs/qc/bin/:${PATH}"

# "Activate" the environment
SHELL ["conda", "run", "--no-capture-output", "-n", "qc", "/bin/bash", "-c"]

# Run
CMD "bash"