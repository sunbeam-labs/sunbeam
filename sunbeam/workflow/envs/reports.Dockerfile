FROM python:3.12-slim

# Setup
WORKDIR /home/reports_env

# Install environment
RUN pip install numpy pandas

RUN echo "Python: $(python --version), Numpy $(pip show numpy | grep 'Version: '), Pandas $(pip show pandas | grep 'Version: ')" > installed_packages.txt

# Run
CMD ["bash"]