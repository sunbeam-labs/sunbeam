#!/bin/bash
set -e
#set -x

PREFIX=$HOME/miniconda3

SUNBEAM_ENV_NAME=${1-sunbeam}
OUTPUT=${2-/dev/stdout}

install_conda () {
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $PREFIX >> $OUTPUT
    export PATH=$PATH:$PREFIX/bin
    command -v conda >/dev/null 2>&1 || { echo "Conda still isn't on the path, try installing manually"; exit 1; }
}

command -v conda >/dev/null 2>&1 || { echo "Conda not installed, installing now"; install_conda; }

conda config --add channels r
conda config --add channels bioconda
conda config --add channels eclarke
conda config --add channels conda-forge

# Don't create the enviroment if it already exists
conda env list | grep -Fxq $SUNBEAM_ENV_NAME || {
    conda create --name=$SUNBEAM_ENV_NAME --file=conda-requirements.txt --yes >> $OUTPUT
    source activate $SUNBEAM_ENV_NAME
    pip install --editable . >> $OUTPUT
    pip install git+https://github.com/eclarke/decontam.git >> $OUTPUT
    echo "Sunbeam successfully installed.";
}

echo "To get started, ensure ${PREFIX}/bin is in your path and run 'source activate $SUNBEAM_ENV_NAME'"

