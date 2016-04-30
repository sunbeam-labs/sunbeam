# Sunbeam: an iridescent HTS pipeline

<img src="http://i.imgur.com/VW3pvQM.jpg" width=320>

Named after the iridescent scales of the sunbeam snake, _Xenopeltis unicolor_, this is a [Snakemake](https://bitbucket.org/snakemake/snakemake/)-based pipeline for high-throughput sequencing. It handles a variety of things that can hopefully enable the reproducible analysis of microbiome (esp. virome) data.

## Installation

Sunbeam uses conda; follow the installation instructions for Miniconda3. Then:

```sh
conda config --add channel r
conda config --add channel biopython
conda env create -f xenopeltis_env.yaml
```
