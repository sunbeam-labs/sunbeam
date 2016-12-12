#!/bin/bash
#set -e
#set -x

# Ensure we can activate the environment
export PATH=$PATH:$HOME/miniconda3/bin

# Activate the sunbeam environment
source activate sunbeam
command -v snakemake

# Temporary
pip install git+https://github.com/zhaoc1/decontam.git

pushd tests

mkdir -p local

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

popd

# Running snakemake
echo "Now testing the snakemake: "
snakemake --configfile=tests/test_config.yml

# Here we just check to ensure it hits the expected genome
echo "Now checking whether we hit the expected genome:"
grep 'NC_006347.1' tests/sunbeam_output/annotation/summary/dummybfragilis.tsv  && return 0 || return 1
