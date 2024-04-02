FROM python:3.12-slim

# Setup
WORKDIR /home/reports_env

# Install environment
RUN pip install numpy pandas

RUN echo "Python: $(python --version)\nNumpy $(pip show numpy | grep 'Version: ')\nPandas $(pip show pandas | grep 'Version: ')" > installed_packages.txt

# Run
CMD ["bash"]