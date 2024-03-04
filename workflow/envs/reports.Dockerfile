FROM python:3.12-slim

# Setup
WORKDIR /home/reports_env

# Install environment
RUN pip install numpy pandas

# Run
CMD ["bash"]