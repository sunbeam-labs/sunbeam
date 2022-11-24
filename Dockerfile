FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="201ccbf210e4ee926c712e6cc1786132a8ffba31966a6a737f07a729fbe44b93"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: ../envs/qc.yml
#   prefix: /conda-envs/fd1dfdfd80cf2f45d2a8d9fc1883f451
#   name: qc
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - biopython
#     - bwa
#     - cutadapt
#     - fastqc
#     - pysam
#     - rust-bio-tools
#     - samtools
#     - trimmomatic
RUN mkdir -p /conda-envs/fd1dfdfd80cf2f45d2a8d9fc1883f451
COPY ../envs/qc.yml /conda-envs/fd1dfdfd80cf2f45d2a8d9fc1883f451/environment.yaml

# Conda environment:
#   source: ../envs/reports.yml
#   prefix: /conda-envs/ffbe99c905f0d2965589357aef4c2358
#   name: reports
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - bioconda::biopython
#     - bioconda::fastqc
#     - conda-forge::pandas
RUN mkdir -p /conda-envs/ffbe99c905f0d2965589357aef4c2358
COPY ../envs/reports.yml /conda-envs/ffbe99c905f0d2965589357aef4c2358/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/fd1dfdfd80cf2f45d2a8d9fc1883f451 --file /conda-envs/fd1dfdfd80cf2f45d2a8d9fc1883f451/environment.yaml && \
    mamba env create --prefix /conda-envs/ffbe99c905f0d2965589357aef4c2358 --file /conda-envs/ffbe99c905f0d2965589357aef4c2358/environment.yaml && \
    mamba clean --all -y
