# Sunbeam: an iridescent HTS pipeline

![xenopeltis unicolor](https://upload.wikimedia.org/wikipedia/commons/thumb/9/97/Young_sunbeam_snake_from_Thailand.JPG/320px-Young_sunbeam_snake_from_Thailand.JPG) 

Named after the iridescent scales of the sunbeam snake, _Xenopeltis unicolor_, this is a [Snakemake](https://bitbucket.org/snakemake/snakemake/)-based pipeline for high-throughput sequencing. It handles a variety of things that can hopefully enable the reproducible analysis of microbiome (esp. virome) data.

## Installation

Sunbeam uses conda; follow the installation instructions for Miniconda3. Then:

```sh
conda config --add channel r
conda config --add channel biopython
conda env create -f xenopeltis_env.yaml
```
