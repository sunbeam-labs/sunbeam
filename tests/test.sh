#!/bin/bash
# Bash flags: Do not commit to repo with these commented out
set -e # Stop on errors
set -x # Echo commands

# Ensure we can activate the environment
export PATH=$PATH:$HOME/miniconda3/bin

# Activate the sunbeam environment
source activate sunbeam
command -v snakemake

# Temporary
pip install git+https://github.com/eclarke/decontam.git

# Set up paths
ROOT=`pwd`
TEMPDIR=`mktemp -d`
mkdir -p $TEMPDIR/data_files

function cleanup {
    # Remove temporary directory if it exists
    # (must be careful with rm -rf and variables)
    [ -z ${TEMPDIR+x} ] || rm -rf "$TEMPDIR"
}

# Calls cleanup when the script exits
trap cleanup EXIT

# Copy data into the temporary directory
cp -r raw $TEMPDIR
cp -r truncated_taxonomy $TEMPDIR
cp -r indexes $TEMPDIR
python generate_dummy_data.py $TEMPDIR

# Create a version of the config file customized for this tempdir
sed "s|TEMPDIR|$TEMPDIR|g" test_config.yml > $TEMPDIR/tmp_config.yml
sed -i "s|HOME|$HOME|g" $TEMPDIR/tmp_config.yml

pushd $TEMPDIR

pwd
ls -ahl
# Generate testing data: data_files
# mkdir -p data_files
# python generate_dummy_data.py

# # Deploy kranken and blast databases
# bash deploy_kraken_db.sh

# bash deploy_blast_db.sh

popd

# Running snakemake
echo "Now testing snakemake: "
#snakemake --configfile=tests/test_config.yml

# Here we just check to ensure it hits the expected genome
echo "Now checking whether we hit the expected genome:"
#grep 'NC_006347.1' tests/sunbeam_output/annotation/summary/dummybfragilis.tsv  && return 0 || return 1
