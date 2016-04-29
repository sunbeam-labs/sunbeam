# Xenopeltis: an iridescent HTS pipeline

![xenopeltis unicolor](https://upload.wikimedia.org/wikipedia/commons/thumb/d/dd/Sunbeam_Snake_%28Xenopeltis_unicolor%29_%287121228691%29.jpg/320px-Sunbeam_Snake_%28Xenopeltis_unicolor%29_%287121228691%29.jpg) 

Named after the iridescent scales of the sunbeam snake, _Xenopeltis unicolor_, this is a [Snakemake](https://bitbucket.org/snakemake/snakemake/)-based pipeline for high-throughput sequencing. It handles a variety of things that can hopefully enable the reproducible analysis of microbiome (esp. virome) data.

## Installation

Xenopeltis uses conda; follow the installation instructions for Miniconda3. Then:

```sh
conda config --add channel r
conda config --add channel biopython
conda env create -f xenopeltis_env.yaml
```
