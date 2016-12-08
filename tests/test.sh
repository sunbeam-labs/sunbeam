#!/bin/bash
set -e
#set -x

# Ensure we can activate the environment
export PATH=$PATH:$HOME/miniconda3/bin

# Activate the sunbeam environment
source activate sunbeam
command -v snakemake

if [ ! local ]; then
	mkdir local
fi

# Generate testing data: data_files
if [ ! -d data_files ]; then
	mkdir data_files
	python generate_dummy_data.py
fi

# Deploy kranken and blast databases
if [ ! -d local/db/bacteria ]; then
	bash deploy_kraken_db.sh
fi

if [ ! -d local/blast ]; then
	bash deploy_blast_db.sh
fi

# Running snakemake
echo "Go to the top folder and run the following: "
echo "snakemake --configfile=tests/test_config.yml"

# Here we just check to ensure it hits the expected genome
echo "To check whether we hit the expected genome, run the following command:"
echo "cat tests/sunbeam_output/annotation/summary/dummybfragilis.tsv | grep 'NC_006347.1' "
