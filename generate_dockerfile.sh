#!/bin/bash

mv extensions/ extensions_moved/
snakemake all --configfile=docker_config.yml --containerize > Dockerfile || true
mv extensions_moved/ extensions