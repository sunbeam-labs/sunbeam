#!/bin/bash

# Allow conda [de]activate
export PATH=$PATH:$HOME/miniconda3/bin
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

# Dump environment contents
conda activate sunbeam
conda list --export
