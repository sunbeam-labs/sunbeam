FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/sunbeam

RUN mkdir -p etc/
COPY etc/* etc/

RUN mkdir -p extensions/
COPY extensions/.placeholder extensions/

RUN mkdir -p src/sunbeamlib/
COPY src/sunbeamlib/* src/sunbeamlib/

COPY tests/ tests/

COPY workflow/ workflow/

COPY environment.yml install.sh MANIFEST.in pyproject.toml pytest.ini README.md ./

# Install sunbeam
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git vim
RUN ./install.sh -e sunbeam -v

ENV PATH="/opt/conda/envs/sunbeam/bin/:${PATH}"
ENV SUNBEAM_DIR="/home/sunbeam"
ENV SUNBEAM_VER="4.4.0"
ENV SUNBEAM_MIN_MEM_MB="8000"
ENV SUNBEAM_MIN_RUNTIME="60"

# Install conda environments
RUN mkdir -p projects
RUN sunbeam init --data_fp tests/data/reads projects/init
RUN sunbeam run --profile projects/init --conda-create-envs-only --mamba
RUN rm -r projects

# "Activate" the environment
SHELL ["conda", "run", "-n", "sunbeam", "/bin/bash", "-c"]

RUN echo "Python: $(python --version)\nSnakemake: $(snakemake --version)\nConda: $(conda --version)" > installed_packages.txt

# Run
CMD "bash"