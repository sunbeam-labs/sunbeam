<img src="docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# Sunbeam: a robust, extensible metagenomic sequencing pipeline 

[![CircleCI](https://circleci.com/gh/sunbeam-labs/sunbeam/tree/dev.svg?style=shield)](https://circleci.com/gh/sunbeam-labs/sunbeam/tree/dev) [![Documentation Status](https://readthedocs.org/projects/sunbeam/badge/?version=latest)](http://sunbeam.readthedocs.io/en/latest/?badge=latest)


Sunbeam is a pipeline written in [snakemake](http://snakemake.readthedocs.io)
that simplifies and automates many of the steps in metagenomic sequencing
analysis. It uses [conda](http://conda.io) to manage dependencies, so it
doesn't have pre-existing dependencies or admin privileges, and can be deployed
on most Linux and Mac workstations and clusters.

Sunbeam currently automates the following tasks:

* Quality control, including adaptor trimming, host read removal, and quality
  filtering;
* Taxonomic assignment of reads to databases using [Kraken](https://github.com/DerrickWood/kraken);
* Assembly of reads into contigs using [Megahit](https://github.com/voutcn/megahit);
* Contig annotation using BLAST[n/p/x];
* Mapping of reads to target genomes; and
* ORF prediction using [Prodigal](https://github.com/hyattpd/Prodigal).

Sunbeam was designed to be modular and extensible. We have a few pre-built
extensions available that handle visualization tasks, including contig
assembly graphs, read alignments, and taxonomic classifications.

To get started, see our [documentation!](https://sunbeam.readthedocs.io)


