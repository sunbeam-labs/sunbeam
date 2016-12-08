#!/bin/bash
set -e
set -x

# Ensure we can activate the environment
export PATH=$PATH:$HOME/miniconda3/bin

# Activate the sunbeam environment
source activate sunbeam
command -v snakemake

mkdir local
mkdir data_files

# Generate testing data: data_files
python generate_dummy_data.py

# Deploy kranken and blast databases
bash deploy_kraken_db.sh
bash deploy_blast_db.sh 


# Running snakemake
echo $"Go to the top folder and run the following: "
echo $"snakemake --configfile=tests/test_config.yml"

# Here we just check to ensure it hits the expected genome
cat tests/sunbeam_output/annotation/summary/dummybfragilis.tsv | grep "NC_006347.1" && return 0 || return 1

