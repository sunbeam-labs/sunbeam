#!/bin/bash

# Allow conda [de]activate
export PATH=$PATH:$HOME/miniconda3/bin
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

# Dump environment contents
TAG=$(git describe --tag)
conda activate sunbeam${TAG:1}
conda list --explicit
