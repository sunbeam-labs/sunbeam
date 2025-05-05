# Pin to version because komplexity is 7 years old and doesn't work with latest Rust
FROM rust:1.80.1-slim

# Setup
WORKDIR /home/komplexity_env

# Install environment
RUN apt-get update && apt-get install -y git
RUN git clone https://github.com/eclarke/komplexity
RUN cd komplexity && cargo install

# Run
CMD ["bash"]