#!/bin/bash
set -e
#set -x

PREFIX=$HOME/miniconda3

install_conda () {
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $PREFIX
    export PATH=$PATH:$PREFIX/bin
    command -v conda >/dev/null 2>&1 || { echo "Conda still isn't on the path, try installing manually"; exit 1; }
}

command -v conda >/dev/null 2>&1 || { echo "Conda not installed, installing now"; install_conda; }

conda config --add channels r
conda config --add channels bioconda
conda config --add channels eclarke

# Don't create the enviroment if it already exists
conda env list | grep -q sunbeam || {
    conda create --name=sunbeam --file=requirements.txt --yes;
    echo "Sunbeam successfully installed.";
}

echo "To get started, ensure ${PREFIX}/bin is in your path and run 'source activate sunbeam'"
