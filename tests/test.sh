#!/bin/bash
set -e
set -x

# Ensure we can activate the environment
export PATH=$PATH:$HOME/miniconda3/bin
source activate sunbeam
command -v snakemake

if [ ! -d tests/sunbeam-data ]; then
    git clone --depth=1 https://github.com/eclarke/sunbeam-data tests/sunbeam-data
fi

mkdir local/db

# Running snakemake without options should give a list of samples
# Here we just check to ensure it outputted a sample name
snakemake --configfile=tests/test-config.yml | grep HUP3D04 && return 0 || return 1

