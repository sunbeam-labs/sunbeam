# Sunbeam: an iridescent virome pipeline

<img src="http://i.imgur.com/VW3pvQM.jpg" width=320>

Named after the iridescent scales of the sunbeam snake, _Xenopeltis unicolor_, this is a [Snakemake](https://bitbucket.org/snakemake/snakemake/)-based pipeline for identifying viral artifacts from high-throughput sequencing reads. 

Right now, the pipeline handles:

- quality control (Trimmomatic, fastqc)
- host read filtering (PMCP _decontam_)
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

- quality control
- read classification
- contig assembly 
- contig annotation

Each section has Snakemake rules defined in the relevant `rules/` subfolder. The
main Snakefile on the root folder loads these rules and assembles the list of
samples to work on.

### Defining samples to work on

- data_fp: the path to the raw files, usually in *.fastq or *.fastq.gz format
- output_fp: the top directory under which the results will be stored
- filename_fmt: how the filenames are stored in `data_fp`
  
    For instance, if you have a list of samples like such:

	- `HUP3D01_L001_R1.fastq.gz`
	- `HUP3D01_L001_R2.fastq.gz`
	- `HUP3D02_L001_R1.fastq.gz`
	- `HUP3D02_L001_R2.fastq.gz`
	- `...`
	
	The sample would be the part defined `HUP3D01`, `HUP3D02`, etc; the segment
	defining the read pair would be `R1` and `R2`, and the intermediate part
	`L001` doesn't change between samples. Thus, the way you'd specify the
	`filename_fmt` parameter would be `{sample}_L001_{rp}.fastq.gz`.
