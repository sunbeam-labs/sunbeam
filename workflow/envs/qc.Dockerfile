FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/qc_env

COPY workflow/envs/qc.yml ./

# Install environment
#SHELL ["conda", "run", "--no-capture-output", "-n", "qc", "/bin/bash", "-c"]
RUN conda env create --file qc.yml --name qc

ENV PATH="${PATH}:/opt/conda/envs/qc/bin/"

# Run
CMD ["conda", "run", "--no-capture-output", "-n", "qc", "/bin/bash", "-c"]