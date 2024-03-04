FROM condaforge/mambaforge:latest

# Setup
RUN useradd --create-home --shell /bin/bash app_user
WORKDIR /home/app_user

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
    apt-get install -y git
RUN ./install.sh -e sunbeam -v

ENV PATH="/opt/conda/envs/sunbeam/bin/:${PATH}"
ENV SUNBEAM_DIR="/home/app_user"
ENV SUNBEAM_VER="4.4.0"
ENV SUNBEAM_MIN_MEM_MB="8000"
ENV SUNBEAM_MIN_RUNTIME="60"

# Install conda environments
#RUN mkdir -p projects
#SHELL ["conda", "run", "--no-capture-output", "-n", "sunbeam_env", "/bin/bash", "-c"]
#RUN sunbeam init projects/init
#RUN sunbeam run --profile projects/init --conda-create-envs-only

# "Activate" the environment
USER app_user
SHELL ["conda", "run", "-n", "sunbeam", "/bin/bash", "-c"]
#RUN bash /opt/conda/envs/sunbeam/etc/conda/activate.d/env_vars.sh

# Run
CMD "bash"