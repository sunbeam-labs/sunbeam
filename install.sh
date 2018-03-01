#!/bin/bash
set -e
#set -x

PREFIX=$HOME/miniconda3

SUNBEAM_ENV_NAME=${1-sunbeam}
OUTPUT=${2-/dev/stdout}

export PATH=$PATH:$PREFIX/bin

install_conda () {
    wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $PREFIX >> $OUTPUT
    export PATH=$PATH:$PREFIX/bin
    command -v conda >/dev/null 2>&1 || { echo "Conda still isn't on the path, try installing manually"; exit 1; }
}

# Install conda if it doesn't show up on the path
command -v conda >/dev/null 2>&1 || { echo "Conda not installed, installing now"; install_conda; }

# Create the environment if it doesn't yet exist
conda env list | cut -f1 -d' ' | grep -Fxq $SUNBEAM_ENV_NAME || {
    
    conda config --add channels r
    conda config --add channels bioconda
    conda config --add channels eclarke
    conda config --add channels conda-forge
    conda create --name=$SUNBEAM_ENV_NAME --file=conda-requirements.txt --quiet --yes >> $OUTPUT
    
}

# Install sunbeam module 
command -v sunbeam_init >/dev/null 2>&1 || {
    source activate $SUNBEAM_ENV_NAME
    pip install --upgrade --editable . >> $OUTPUT
    command -v conda >/dev/null 2>&1 || { echo "Couldn't install sunbeam; please report this as a bug."; exit 1; }
    echo "Sunbeam successfully installed.";
}

echo "To get started, ensure ${PREFIX}/bin is in your path and run 'source activate $SUNBEAM_ENV_NAME'"

