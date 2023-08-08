#!/bin/bash

# Allow conda [de]activate
export PATH=$PATH:$HOME/miniconda3/bin
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

TAG=$(git describe --tag)
conda activate sunbeam${TAG:1}

DIR=$SUNBEAM_DIR/.snakemake/

ls -la $DIR
tail -n +1 $DIR/*.pin.txt
#for fp in $DIR/*.yaml; do
#    cat $fp
#    ENV=${fp%.yaml}
#    conda list -p $ENV
#done