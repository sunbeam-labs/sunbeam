FROM condaforge/mambaforge:latest

WORKDIR /home/sunbeam

# Install Sunbeam
COPY pyproject.toml .
COPY sunbeam/ sunbeam/
RUN pip install .[dev]

# Pre-create conda environments (caches them in image)
RUN mkdir -p projects && \
    sunbeam init --data_fp tests/data/reads projects/init && \
    sunbeam run --profile projects/init --conda-create-envs-only --mamba && \
    rm -rf projects/

# Set labels
LABEL org.opencontainers.image.title="Sunbeam" \
      org.opencontainers.image.source="https://github.com/sunbeam-labs/sunbeam"

# Set entry point
CMD ["sunbeam"]
