# Sunbeam: an iridescent high-throughput sequencing pipeline 
<img src="http://i.imgur.com/VW3pvQM.jpg" width=240> 

[![Build Status](https://travis-ci.org/eclarke/sunbeam.svg?branch=master)](https://travis-ci.org/eclarke/sunbeam) 
[![Code Health](https://landscape.io/github/eclarke/sunbeam/master/landscape.svg?style=flat)](https://landscape.io/github/eclarke/sunbeam/master)
[![Snakemake](https://img.shields.io/badge/snakemake-≥3.5.2-brightgreen.svg?style=flat)](http://snakemake.bitbucket.org)


Named after the iridescent scales of the sunbeam snake, _Xenopeltis unicolor_,
this is a [Snakemake](https://bitbucket.org/snakemake/snakemake/)-based pipeline
for creating various products from high-throughput sequencing reads.

Right now, the pipeline handles:

- quality control (Trimmomatic, fastqc)
- host read filtering (PMCP _decontam_)
- kmer classification (kraken)
- contig assembly (idba_ud + CAP3)
- nucleotide search and alignment (blastn, bowtie2)
- ORF prediction (MetaGene Annotator)
- protein blast on predicted genes (blastp, blastx)
- summarizing blast results

## Installation

Sunbeam uses [conda](http://conda.pydata.org/miniconda.html). To deploy
automatically, run `install.sh`. The script will download and install Miniconda3
to the home folder automatically if it does not exist, and builds an environment
with the needed dependencies. If you're having issues, make sure to add 
`miniconda3/bin` to your `PATH`.

## High-level overview

The pipeline is separated into sections: 

- quality control (`qc`)
- read classification (`classify`)
- read mapping and alignment visualization (`mapping`)
- contig assembly (`assembly`)
- contig annotation (`annotation`)

Each section has Snakemake rules defined in the relevant `rules/` subfolder. The
main Snakefile on the root folder loads these rules and assembles the list of
samples to work on.

## Using sunbeam

Before you get started, you just need to have a folder in mind that will be your
'project folder'. This should contain a subfolder with your sequence data
in it.

To get started, activate the sunbeam environment via `source activate sunbeam` 
and then run the following command, with the name of your folder instead of
`{path_to_project_folder}` and a custom name instead of `{project_config.yml}`.

```
sunbeam_init {path_to_project_folder} > {project_config.yml}
```

For example, if I have a project about ocean viruses, I would run this:

`sunbeam_init /home/erik/OceanVirome > ocean_virome_config.yml`

Next, edit the resulting config file to point to the various paths for your data and
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

### Defining genomes to align to

The mapping rules will attempt to align all reads from all samples to all
available genome sequences, and will also create some summary statistics and a
visualization of alignments for each genome using [IGV].

At a minimum the mapping configuration requires a path to a directory of genome
files in fasta format (for alignment of reads), and the path to an IGV
executable (for visualizing alignments against each genome).  IGV version
2.3.68 is known to work.  [Xvfb] also must be installed for IGV visualization.

An example for the two required mapping configuration options, if you have a
custom IGV install in `sunbeam/local/IGV`:

    mapping:
      genomes_fp: "/home/erik/OceanVirome/genomes"
      igv_fp: "local/IGV/igv.sh"

Some other settings specific to the mapping rules:

- `keep_unaligned`: A boolean defining if unaligned reads should be kept in
  bowtie2's output, along with aligned reads. (Defaults to False.)
- `igv_prefs`: A dictionary of IGV preferences to apply.  (Defaults to a few
  basic rendering options.)

### Databases

Currently, a number of the annotation and classification rules require specific
copies of their databases to be accessible. If you are not ready to use them,
you can leave their paths blank and the automated checker will not complain
about them being invalid. 

#### Setting up custom Kraken databases

[Kraken](https://ccb.jhu.edu/software/kraken/) is the default read annotation software for Sunbeam.
It is installed as a dependency, so manual installation is not required.
Unfortunately it hasn't been updated in a while and the scripts to build custom databases no longer work with NCBI.

We recommend installing this [fork](https://github.com/taltman/kraken), which updates the database building scripts to use the correct NCBI URLs. 
As a bonus, it also uses `dust` to mask problematic low-complexity regions of the genomes being used to build the database. 
This step is essential to perform by hand if not using this fork; see this page on instructions: https://groups.google.com/forum/#!topic/kraken-users/jjRe21-qyvw

### Running sunbeam

```sh
snakemake --configfile=<my-config.yml> [optional rule]
```

Running sunbeam without specifying a rule will perform, by default:

- quality control
- host read filtering
- read alignment to genome sequences
- read-level classification from given Kraken database
- contig assembly
- gene detection
- blasting of genes and contigs against given databases

Each step of this can be run piecemeal by using the rules starting with "all", such as:

- all_qc: quality control on all reads
- all_decontam: remove all host reads
- all_mapping: align all reads to all genome sequences
- all_classify: classify all reads
- all_assembly: build contigs
- all_annotate: annotate all contigs

Example: decontaminate all host reads with ```snakemake --configfile=my_config.yml all_decontam```

You can also see what samples Sunbeam detected by running ```snakemake --configfile=my_config.yml samples```.


## Updating Sunbeam

To update sunbeam when a new version is released, it's easiest to simply remove the sunbeam conda environment and reinstall:

```shell
conda env remove sunbeam
cd path/to/sunbeam
git pull
bash install.sh
```

## Troubleshooting

- **Q**: When I type `source activate sunbeam`, I get an error.
- **A**: Make sure that miniconda3 is in your `PATH`: type `echo $PATH` and look for miniconda there. If it does not exist, edit your `.profile` or `.bashrc` file in your home directory and add a line like this: 
    
    ```
    export PATH=$PATH:$HOME/miniconda3/bin
    ```

- **Q**: I get an error relating to missing files when running on a cluster, and nothing continues after that point.
- **A**: Increase the waiting time for files to appear by adding `-w 90` to the snakemake command. 

[IGV]: https://software.broadinstitute.org/software/igv
[Xvfb]: https://www.x.org/archive/X11R7.7/doc/man/man1/Xvfb.1.xhtml
