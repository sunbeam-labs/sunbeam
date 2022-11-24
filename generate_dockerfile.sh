#!/bin/bash

mv extensions/ extensions_moved/
snakemake --configfile=$1 --containerize > Dockerfile || true
mv extensions_moved/ extensions
docker build - < Dockerfile