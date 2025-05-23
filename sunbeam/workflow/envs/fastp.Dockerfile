FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/fastp_env

COPY sunbeam/workflow/envs/fastp.yml ./

# Install environment
RUN mamba env create --file fastp.yml --name fastp

ENV PATH="/opt/conda/envs/fastp/bin/:${PATH}"

# "Activate" the environment
SHELL ["conda", "run", "--no-capture-output", "-n", "fastp", "/bin/bash", "-c"]

# Run
CMD "bash"