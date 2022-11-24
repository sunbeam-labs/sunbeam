#!/bin/bash

mv extensions/ extensions_moved/
snakemake all --configfile=projects/WGS-test/sunbeam_config.yml --containerize > Dockerfile || true
mv extensions_moved/ extensions