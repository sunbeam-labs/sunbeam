FROM condaforge/mambaforge:latest

WORKDIR /home/sunbeam

# Install Sunbeam
COPY pyproject.toml .
COPY sunbeam/ sunbeam/
RUN pip install .

LABEL org.opencontainers.image.title="Sunbeam Slim" \
      org.opencontainers.image.source="https://github.com/sunbeam-labs/sunbeam"

# Set entry point
CMD ["sunbeam"]
