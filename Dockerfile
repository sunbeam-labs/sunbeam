FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="a0adbd8232d5365eb11a5b01b898e0f31c6fc417d72fe95e856d5800f23be10d"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: ../../../envs/assembly.yml
#   prefix: /conda-envs/418b67981e43f64151d2840602a45c78
#   name: assembly
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - megahit # https://github.com/sunbeam-labs/sunbeam/issues/213
RUN mkdir -p /conda-envs/418b67981e43f64151d2840602a45c78
COPY ../../../envs/assembly.yml /conda-envs/418b67981e43f64151d2840602a45c78/environment.yaml

# Conda environment:
#   source: ../../../envs/komplexity.yml
#   prefix: /conda-envs/f2039bdd3dcc9f605689b977f5422cfc
#   name: komplexity
#   channels:
#     - eclarke
#     - conda-forge # Required by komplexity recipe
#   dependencies:
#     - komplexity
RUN mkdir -p /conda-envs/f2039bdd3dcc9f605689b977f5422cfc
COPY ../../../envs/komplexity.yml /conda-envs/f2039bdd3dcc9f605689b977f5422cfc/environment.yaml

# Conda environment:
#   source: ../../../envs/qc.yml
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
COPY ../../../envs/qc.yml /conda-envs/fd1dfdfd80cf2f45d2a8d9fc1883f451/environment.yaml

# Conda environment:
#   source: ../../../envs/reports.yml
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
COPY ../../../envs/reports.yml /conda-envs/ffbe99c905f0d2965589357aef4c2358/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/418b67981e43f64151d2840602a45c78 --file /conda-envs/418b67981e43f64151d2840602a45c78/environment.yaml && \
    mamba env create --prefix /conda-envs/f2039bdd3dcc9f605689b977f5422cfc --file /conda-envs/f2039bdd3dcc9f605689b977f5422cfc/environment.yaml && \
    mamba env create --prefix /conda-envs/fd1dfdfd80cf2f45d2a8d9fc1883f451 --file /conda-envs/fd1dfdfd80cf2f45d2a8d9fc1883f451/environment.yaml && \
    mamba env create --prefix /conda-envs/ffbe99c905f0d2965589357aef4c2358 --file /conda-envs/ffbe99c905f0d2965589357aef4c2358/environment.yaml && \
    mamba clean --all -y
