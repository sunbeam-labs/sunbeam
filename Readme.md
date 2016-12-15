# Sunbeam: an iridescent virome pipeline 
<img src="http://i.imgur.com/VW3pvQM.jpg" width=240> 

[![Build Status](https://travis-ci.org/eclarke/sunbeam.svg?branch=master)](https://travis-ci.org/eclarke/sunbeam) 
[![Code Health](https://landscape.io/github/eclarke/sunbeam/master/landscape.svg?style=flat)](https://landscape.io/github/eclarke/sunbeam/master)
[![Snakemake](https://img.shields.io/badge/snakemake-≥3.5.2-brightgreen.svg?style=flat)](http://snakemake.bitbucket.org)


Named after the iridescent scales of the sunbeam snake, _Xenopeltis unicolor_,
this is a [Snakemake](https://bitbucket.org/snakemake/snakemake/)-based pipeline
for identifying viral artifacts from high-throughput sequencing reads.

Right now, the pipeline handles:

- quality control (Trimmomatic, fastqc)
- host read filtering (PMCP _decontam_)
- kmer classification (kraken, CLARK)
- visualization (krona)
- contig assembly (idba_ud, minimo)
- nucleotide search and alignment (blastn, bowtie2)
- ORF prediction (MetaGene Annotator)
- protein blast on predicted genes (blastp, blastx)
- summarizing blast results

## Installation

Sunbeam uses [conda](http://conda.pydata.org/miniconda.html). To deploy
automatically, run `install.sh`. The script will download and install Miniconda3
to the home folder automatically if it does not exist, and builds an environment
with the needed dependencies.

## High-level overview

The pipeline is separated into sections: 

- quality control (`qc`)
- read classification (`classify`)
- contig assembly (`assembly`)
- contig annotation (`annotation`)

Each section has Snakemake rules defined in the relevant `rules/` subfolder. The
main Snakefile on the root folder loads these rules and assembles the list of
samples to work on.

## Using sunbeam

To get started, make a copy of the `tests/test-config.yml` file in the same
directory as the `Snakefile` and name it something like `my-config.yml`.

Next, edit the config file to point to the various paths for your data and
output directories, as shown below:

### Defining samples to work on

- `data_fp`: the path to the raw files, usually in *.fastq or *.fastq.gz format
- `output_fp`: the top directory under which the results will be stored
- `filename_fmt`: how the filenames are stored in `data_fp`
  
For instance, if you have a list of samples like such:

- `HUP3D01_L001_R1.fastq.gz`
- `HUP3D01_L001_R2.fastq.gz`
- `HUP3D02_L001_R1.fastq.gz`
- `HUP3D02_L001_R2.fastq.gz`
- `...`
	
The sample would be the part defined `HUP3D01`, `HUP3D02`, etc; the segment
defining the read pair would be `R1` and `R2`, and the intermediate part `L001`
doesn't change between samples. Thus, the way you'd specify the `filename_fmt`
parameter would be `{sample}_L001_{rp}.fastq.gz`.

### Databases

Currently, a number of the annotation and classification rules require specific
copies of their databases to be accessible. If you are not ready to use them,
you can leave their paths as blank and the automated checker will not complain
about them being invalid. 

### Running sunbeam

```sh
snakemake --configfile=<my-config.yml>
```

will print out a list of samples that it found according to your specified
config file and invites you to peruse the list of rules by executing 
`snakemake --configfile=<my-config.yml> --list`.
