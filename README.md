<img src="docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# Sunbeam: a robust, extensible metagenomic sequencing pipeline 

[![CircleCI](https://circleci.com/gh/sunbeam-labs/sunbeam/tree/main.svg?style=shield)](https://circleci.com/gh/sunbeam-labs/sunbeam/tree/main) [![Super-Linter](https://github.com/sunbeam-labs/sunbeam/actions/workflows/linter.yml/badge.svg)](https://github.com/sunbeam-labs/sunbeam/actions/workflows/linter.yml) [![Conda Envs Status](https://byob.yarr.is/sunbeam-labs/sunbeam/env_check)] [![Documentation Status](https://readthedocs.org/projects/sunbeam/badge/?version=stable)](https://sunbeam.readthedocs.io/en/stable/?badge=stable) [![DOI:10.1186/s40168-019-0658-x](https://img.shields.io/badge/Published%20in-Microbiome-1abc9c.svg)](https://doi.org/10.1186/s40168-019-0658-x)

Sunbeam is a pipeline written in [snakemake](http://snakemake.readthedocs.io)
that simplifies and automates many of the steps in metagenomic sequencing
analysis. It uses [conda](http://conda.io) to manage dependencies, so it
doesn't have pre-existing dependencies or admin privileges, and can be deployed
on most Linux workstations and clusters. To read more, check out [our paper
in Microbiome](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-019-0658-x).

Sunbeam currently automates the following tasks:

* Quality control, including adaptor trimming, host read removal, and quality
  filtering;
* Taxonomic assignment of reads to databases using [Kraken](https://github.com/DerrickWood/kraken);
* Assembly of reads into contigs using [Megahit](https://github.com/voutcn/megahit);
* Contig annotation using BLAST[n/p/x];
* Mapping of reads to target genomes; and
* ORF prediction using [Prodigal](https://github.com/hyattpd/Prodigal).

Sunbeam was designed to be modular and extensible. Some extensions have been built for:

- [IGV](https://github.com/sunbeam-labs/sbx_igv) for viewing read alignments
- [KrakenHLL](https://github.com/zhaoc1/sbx_krakenhll), an alternate read classifier
- [Kaiju](https://github.com/sunbeam-labs/sbx_kaiju), a read classifier that uses BWA rather than kmers
- [Anvi'o](https://github.com/sunbeam-labs/sbx_anvio), a downstream analysis pipeline that does lots of stuff!

More extensions can be found at the extension page: https://github.com/sunbeam-labs.

**To get started, see our [documentation](http://sunbeam.readthedocs.io)!**

If you use the Sunbeam pipeline in your research, please cite: 

EL Clarke, LJ Taylor, C Zhao, A Connell, J Lee, FD Bushman, K Bittinger. Sunbeam: an
extensible pipeline for analyzing metagenomic sequencing experiments. *Microbiome* 7:46 (2019)

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

