#!/bin/bash
# $1 is the version for this release (e.g. 4.0.0 or 4.1.3-rc.2)

snakemake --configfile dev_scripts/sunbeam_config.yml --archive tmp.tar.gz
rm tmp.tar.gz
tar -czvf sunbeam$1.tar.gz .snakemake/conda-archive/ etc/ extensions/.placeholder sunbeamlib/ tests/ workflow/ environment.yml install.sh Readme.md setup.py
