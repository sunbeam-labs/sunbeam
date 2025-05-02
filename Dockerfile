FROM condaforge/mambaforge:latest

WORKDIR /home/sunbeam

# Install Sunbeam
COPY pyproject.toml .
COPY sunbeam/ sunbeam/
RUN pip install .

# Pre-create conda environments (caches them in image)
COPY tests/data/reads/ reads/
RUN mkdir -p projects && \
    sunbeam init --data_fp reads projects/init && \
    sunbeam run --profile projects/init --conda-create-envs-only --mamba && \
    rm -r projects reads

LABEL org.opencontainers.image.title="Sunbeam" \
      org.opencontainers.image.source="https://github.com/sunbeam-labs/sunbeam"

# Set entry point
CMD ["sunbeam"]
