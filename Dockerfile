FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="12452d6af82ca1066b553b1905c85b62405942009ebafb5e897d30b84646e799"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/cutadapt.yml
#   prefix: /conda-envs/08b0272f8c744b8bb162030774cf9917
#   channels:
#     - bioconda
#   dependencies:
#     - cutadapt
#     #- python =3.12.0
#   name: cutadapt
RUN mkdir -p /conda-envs/08b0272f8c744b8bb162030774cf9917
COPY workflow/envs/cutadapt.yml /conda-envs/08b0272f8c744b8bb162030774cf9917/environment.yaml

# Conda environment:
#   source: workflow/envs/komplexity.yml
#   prefix: /conda-envs/cc72134ebb57f052a6ac7d5a084f9bf6
#   channels:
#     - eclarke
#     - conda-forge
#   dependencies:
#     - komplexity
#     - python =3.12.0
#   name: komplexity
RUN mkdir -p /conda-envs/cc72134ebb57f052a6ac7d5a084f9bf6
COPY workflow/envs/komplexity.yml /conda-envs/cc72134ebb57f052a6ac7d5a084f9bf6/environment.yaml

# Conda environment:
#   source: workflow/envs/qc.yml
#   prefix: /conda-envs/2e2d3aeac664927deed435e9fde5bf30
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - bwa
#     - fastqc
#     - trimmomatic
#     - python =3.12.0
#   name: qc
RUN mkdir -p /conda-envs/2e2d3aeac664927deed435e9fde5bf30
COPY workflow/envs/qc.yml /conda-envs/2e2d3aeac664927deed435e9fde5bf30/environment.yaml

# Conda environment:
#   source: workflow/envs/reports.yml
#   prefix: /conda-envs/60172346a36951e9f4497f0c766b8d74
#   channels:
#     - conda-forge
#   dependencies:
#     - numpy
#     - pandas
#     - python =3.12.0
#   name: reports
RUN mkdir -p /conda-envs/60172346a36951e9f4497f0c766b8d74
COPY workflow/envs/reports.yml /conda-envs/60172346a36951e9f4497f0c766b8d74/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/08b0272f8c744b8bb162030774cf9917 --file /conda-envs/08b0272f8c744b8bb162030774cf9917/environment.yaml && \
    mamba env create --prefix /conda-envs/cc72134ebb57f052a6ac7d5a084f9bf6 --file /conda-envs/cc72134ebb57f052a6ac7d5a084f9bf6/environment.yaml && \
    mamba env create --prefix /conda-envs/2e2d3aeac664927deed435e9fde5bf30 --file /conda-envs/2e2d3aeac664927deed435e9fde5bf30/environment.yaml && \
    mamba env create --prefix /conda-envs/60172346a36951e9f4497f0c766b8d74 --file /conda-envs/60172346a36951e9f4497f0c766b8d74/environment.yaml && \
    mamba clean --all -y
