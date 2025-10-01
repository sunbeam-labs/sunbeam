FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/decontam_env

COPY sunbeam/workflow/envs/decontam.yml ./

# Install environment
RUN mamba env create --file decontam.yml --name decontam

ENV PATH="/opt/conda/envs/decontam/bin/:${PATH}"

# "Activate" the environment
SHELL ["conda", "run", "--no-capture-output", "-n", "decontam", "/bin/bash", "-c"]

# Run
CMD "bash"