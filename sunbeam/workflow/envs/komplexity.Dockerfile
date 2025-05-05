FROM rust:slim

# Setup
WORKDIR /home/komplexity_env

# Install environment
RUN apt-get update && apt-get install -y git
RUN git clone https://github.com/eclarke/komplexity
RUN cd komplexity && cargo install

RUN echo "Rustc: $(rustc --version), Komplexity $(kz --version)" > installed_packages.txt

# Run
CMD ["bash"]