FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="7699e86cfbb77f2215ae4988c1451a0b8f80d8be28a977bb92c80bde4ed485e9"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/annotation.yml
#   prefix: /conda-envs/8d6dce10cf0795dd3b072d442d840b64
#   name: annotation
#   channels:
#     - bioconda
#   dependencies:
#     - prodigal
RUN mkdir -p /conda-envs/8d6dce10cf0795dd3b072d442d840b64
COPY workflow/envs/annotation.yml /conda-envs/8d6dce10cf0795dd3b072d442d840b64/environment.yaml

# Conda environment:
#   source: workflow/envs/assembly.yml
#   prefix: /conda-envs/418b67981e43f64151d2840602a45c78
#   name: assembly
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - megahit # https://github.com/sunbeam-labs/sunbeam/issues/213
RUN mkdir -p /conda-envs/418b67981e43f64151d2840602a45c78
COPY workflow/envs/assembly.yml /conda-envs/418b67981e43f64151d2840602a45c78/environment.yaml

# Conda environment:
#   source: workflow/envs/komplexity.yml
#   prefix: /conda-envs/f2039bdd3dcc9f605689b977f5422cfc
#   name: komplexity
#   channels:
#     - eclarke
#     - conda-forge # Required by komplexity recipe
#   dependencies:
#     - komplexity
RUN mkdir -p /conda-envs/f2039bdd3dcc9f605689b977f5422cfc
COPY workflow/envs/komplexity.yml /conda-envs/f2039bdd3dcc9f605689b977f5422cfc/environment.yaml

# Conda environment:
#   source: workflow/envs/qc.yml
#   prefix: /conda-envs/9b43e0e6246759407dc66c30b41c518d
#   name: qc
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - conda-forge::biopython
#     - bwa
#     - cutadapt
#     - fastqc
#     - pysam
#     - rust-bio-tools
#     - samtools
#     - trimmomatic
#     - python=3.10
RUN mkdir -p /conda-envs/9b43e0e6246759407dc66c30b41c518d
COPY workflow/envs/qc.yml /conda-envs/9b43e0e6246759407dc66c30b41c518d/environment.yaml

# Conda environment:
#   source: workflow/envs/reports.yml
#   prefix: /conda-envs/bb8af9913d098aa435b66c7fd30d2ef3
#   name: reports
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - conda-forge::biopython
#     - fastqc
#     - pandas
#     - python=3.10
RUN mkdir -p /conda-envs/bb8af9913d098aa435b66c7fd30d2ef3
COPY workflow/envs/reports.yml /conda-envs/bb8af9913d098aa435b66c7fd30d2ef3/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/8d6dce10cf0795dd3b072d442d840b64 --file /conda-envs/8d6dce10cf0795dd3b072d442d840b64/environment.yaml && \
    mamba env create --prefix /conda-envs/418b67981e43f64151d2840602a45c78 --file /conda-envs/418b67981e43f64151d2840602a45c78/environment.yaml && \
    mamba env create --prefix /conda-envs/f2039bdd3dcc9f605689b977f5422cfc --file /conda-envs/f2039bdd3dcc9f605689b977f5422cfc/environment.yaml && \
    mamba env create --prefix /conda-envs/9b43e0e6246759407dc66c30b41c518d --file /conda-envs/9b43e0e6246759407dc66c30b41c518d/environment.yaml && \
    mamba env create --prefix /conda-envs/bb8af9913d098aa435b66c7fd30d2ef3 --file /conda-envs/bb8af9913d098aa435b66c7fd30d2ef3/environment.yaml && \
    mamba clean --all -y
