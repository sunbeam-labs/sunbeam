<img src="https://github.com/sunbeam-labs/sunbeam/blob/main/docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# Sunbeam: a robust, extensible metagenomic sequencing pipeline 

[![Tests](https://github.com/sunbeam-labs/sunbeam/actions/workflows/tests.yml/badge.svg)](https://github.com/sunbeam-labs/sunbeam/actions/workflows/tests.yml)
[![Documentation Status](https://readthedocs.org/projects/sunbeam/badge/?version=stable)](https://sunbeam.readthedocs.io/en/stable/?badge=stable)
[![PyPI](https://badge.fury.io/py/sunbeamlib.svg)](https://pypi.org/project/sunbeamlib/)
[![Bioconda](https://anaconda.org/bioconda/sunbeamlib/badges/downloads.svg)](https://anaconda.org/bioconda/sunbeamlib/)
[![DockerHub](https://img.shields.io/docker/pulls/sunbeamlabs/sunbeam)](https://hub.docker.com/repository/docker/sunbeamlabs/sunbeam/)
[![DOI:10.1186/s40168-019-0658-x](https://img.shields.io/badge/Published%20in-Microbiome-1abc9c.svg)](https://doi.org/10.1186/s40168-019-0658-x)

Sunbeam is a pipeline written in [snakemake](http://snakemake.readthedocs.io) that simplifies and automates many of the steps in metagenomic sequencing analysis. It can be installed through PyPi (pip), Bioconda (conda), DockerHub (docker), or GitHub (git). Sunbeam was designed to be modular and extensible, allowing anyone to build off the core functionality. To read more, check out [our paper in Microbiome](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-019-0658-x).

Sunbeam currently automates the following tasks:

* Quality control, including adapter trimming, host read removal, and quality filtering;
* Taxonomic assignment of reads to databases using [Kraken](https://github.com/DerrickWood/kraken) ([sbx_kraken](https://github.com/sunbeam-labs/sbx_kraken));
* Assembly of reads into contigs using [Megahit](https://github.com/voutcn/megahit) ([sbx_assembly](https://github.com/sunbeam-labs/sbx_assembly));
* Contig annotation using BLAST[n/p/x] and Diamond ([sbx_assembly](https://github.com/sunbeam-labs/sbx_assembly));
* Mapping to reference genomes ([sbx_mapping](https://github.com/sunbeam-labs/sbx_mapping))
* ORF prediction using [Prodigal](https://github.com/hyattpd/Prodigal) ([sbx_assembly](https://github.com/sunbeam-labs/sbx_assembly)).

More extensions can be found at https://github.com/sunbeam-labs.

**To get started, see our [documentation](http://sunbeam.readthedocs.io)!**

If you use the Sunbeam pipeline in your research, please cite: 

EL Clarke, LJ Taylor, C Zhao, A Connell, J Lee, FD Bushman, K Bittinger. Sunbeam: an extensible pipeline for analyzing metagenomic sequencing experiments. *Microbiome* 7:46 (2019)

See how people are using Sunbeam:

- Shi, Z *et al.* [Segmented Filamentous Bacteria Prevent and Cure Rotavirus Infection.](https://www.sciencedirect.com/science/article/pii/S0092867419310797) *Cell* 179, 644-658.e13 (2019).
- Abbas, AA *et al.* [Redondoviridae, a Family of Small, Circular DNA Viruses of the Human Oro-Respiratory Tract Associated with Periodontitis and Critical Illness.](https://www.sciencedirect.com/science/article/pii/S1931312819301714) *Cell Host Microbe* 25, 719–729 (2019).
- Leiby, JS *et al.* [Lack of detection of a human placenta microbiome in samples from preterm and term deliveries.](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0575-4) *Microbiome* 6, 1–11 (2018).

------

### Contributors

- Erik Clarke ([@eclarke](https://github.com/eclarke))
- Chunyu Zhao ([@zhaoc1](https://github.com/zhaoc1))
- Jesse Connell ([@ressy](https://github.com/ressy))
- Louis Taylor ([@louiejtaylor](https://github.com/louiejtaylor))
- Charlie Bushman ([@ulthran](https://github.com/ulthran))
- Kyle Bittinger ([@kylebittinger](https://github.com/kylebittinger))

