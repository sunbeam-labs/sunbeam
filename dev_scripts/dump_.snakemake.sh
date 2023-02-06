#!/bin/bash

DIR=$SUNBEAM_DIR/.snakemake/

ls -la $DIR
for fp in $DIR/*.yaml; do
    cat $fp
    ENV=${fp%.yaml}
    conda list -p $ENV
done