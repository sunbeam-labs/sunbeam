FROM python:3.12-slim

# Setup
WORKDIR /home/cutadapt_env

# Install environment
RUN pip install cutadapt

RUN echo "Python: $(python --version), Cutadapt $(pip show cutadapt | grep 'Version: ')" > installed_packages.txt

# Run
CMD ["bash"]