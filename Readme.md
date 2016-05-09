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

Sunbeam uses conda; follow the installation instructions for Miniconda3. Then:

```sh
conda config --add channel r
conda config --add channel biopython
conda env create -f xenopeltis_env.yaml
```

### External dependencies

The growing list of external dependencies that need to be installed
because they are unavailable through Conda:

- Trimmomatic
- idba_ud, modified to handle longer reads
- minimo
- MetaGene Annotator (mga)

## High-level overview

The pipeline is separated into three sections: quality control (qc), read
assembly, and read/contig annotation. Each section has Snakemake rules defined
in the relevant folder, along with section-specific config files. The main
Snakefile on the root folder loads these rules and assembles the list of samples
to work on.

### Defining samples to work on

- data_fp: the path to the raw files, usually in *.fastq or *.fastq.gz format
- output_fp: the 
