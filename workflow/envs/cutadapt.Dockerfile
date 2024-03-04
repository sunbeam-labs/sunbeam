FROM python:3.12-slim

# Setup
WORKDIR /home/cutadapt_env

# Install environment
RUN pip install cutadapt

# Run
CMD ["bash"]