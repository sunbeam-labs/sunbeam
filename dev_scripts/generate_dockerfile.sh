#!/bin/bash

snakemake all --configfile=dev_scripts/sunbeam_config.yml --containerize > Dockerfile || true