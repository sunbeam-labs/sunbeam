FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/qc_env

COPY workflow/envs/qc.yml ./

# Install environment
RUN mamba env create --file qc.yml --name qc

ENV PATH="/opt/conda/envs/qc/bin/:${PATH}"

# "Activate" the environment
SHELL ["conda", "run", "--no-capture-output", "-n", "qc", "/bin/bash", "-c"]

RUN echo "Python: $(python --version), Conda: $(conda --version), BWA: $(bwa 2>&1 >/dev/null | grep 'Version: '), FastQC: $(fastqc --version), Trimmomatic: $(trimmomatic -version)" > installed_packages.txt

# Run
CMD "bash"