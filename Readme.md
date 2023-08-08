<img src="docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# Sunbeam: a robust, extensible metagenomic sequencing pipeline 

[![CircleCI](https://circleci.com/gh/sunbeam-labs/sunbeam/tree/dev.svg?style=shield)](https://circleci.com/gh/sunbeam-labs/sunbeam/tree/dev) [![Super-Linter](https://github.com/sunbeam-labs/sunbeam/actions/workflows/linter.yml/badge.svg)](https://github.com/sunbeam-labs/sunbeam/actions/workflows/linter.yml) [![Documentation Status](https://readthedocs.org/projects/sunbeam/badge/?version=stable)](https://sunbeam.readthedocs.io/en/stable/?badge=stable) [![DOI:10.1186/s40168-019-0658-x](https://img.shields.io/badge/Published%20in-Microbiome-1abc9c.svg)](https://doi.org/10.1186/s40168-019-0658-x)

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

### Changelog:

#### v4.0.0 (July 29, 2023)

- Many components of sunbeam have been moved extensions to minimize installation and usage issues with the core pipeline, qc and decontam are the remaining core components
- Deprecated --cluster-config snakemake option replaced with profiles, also reduces sunbeam's command line interface by moving options to the profile
- Improved logging and reporting standards for sunbeam and extensions
- Restructured to snakemake standards for tarball installation methods
- Switch test suite from bash to pytest to increase ease and efficiency of testing
- Updated documentation
- Improved filter_reads efficiency on large hostid sets and catch errors in sequence header matching

#### v3.1.1 (October 19, 2022)

 - Upgrade manage-version.sh to work with any version (before or after 3.1) and always keep itself on the latest version
 - Allow multiple targets with `sunbeam run --target_list [TARGET, ...]`

#### v3.1.0 (October 14, 2022)

 - Upgrade snakemake to v7.15.1, this fixed two issues with running jobs on the cluster
 - Add manage-versions.sh which manages versions of sunbeam automatically for users in a single repo
 - Fix build errors for previous versions of docs introduced when trying to add versioning info from sunbeamlib
 - Improve the install script's decision making on when to use mamba vs conda and add informative errors
 - Remove 'str' calls from snakemake that are no longer necessary (snakemake <5.7 requirement)
 - Add 'Sunbeam Commands', 'Software Structure', 'install.sh', and 'manage-versions.sh' pages to docs

#### v3.0.1 (August 11, 2022)

 - Add biom-format to Kraken env
 - Fix sunbeam init not picking up extension configs
 - Fix multiqc naming issue that seemed to arise across different multiqc versions

#### v3.0.0 (June 27, 2022)

 - Support use of .smk file extensions in Sunbeam extensions (in addition to .rules)
 - Making use of snakemake's builtin features for environment management to separate dependencies and shrink environments
 - Support [mamba](https://github.com/mamba-org/mamba) as an alternate package dependency solver at install time, for faster installs
 - New command `sunbeam extend` to automatically install Sunbeam extensions! Use like `sunbeam extend https://github.com/sunbeam-labs/sbx_report`
 - `sunbeam init` and `sunbeam config update` now add options for extensions you've installed to your default config file! (#247)
 - Updated the path to the Illumina adapter sequences from hardcoded to templated (fixes #150 and #152)
 - Use the updated kraken2 classifier instead of kraken
 - Update other dependencies (trimmomatic -> 0.3.9; grabseqs -> 0.6.1; snakemake -> <5.7.0)
 - Use diamond instead of blastx/p for a significant speed increase

#### v2.1.0 (November 26, 2019)

 - Added a build manifest, which is run every time on integration testing and can be fed into conda by users to install the most recent successful dependencies
 - Updates to documentation (#169, #230, #231)
 - Fix missing samtools (#224)
 - Integration test updates to schedule weekly builds (#222)
 - Fix issues with old paired-end illumina adapters (#221)
 - Script updates to use conda commands instead of source commands (#220)
 - Add h5py package explicitly to avoid dependency metadata problem (#219)
 - Add multiQC to build QC report (#203)
 - Use multithreading for cutadapt in QC (#202)
 - Correct conda channel priority during install (#201)
 - Update documentation to spell out requirements (#199)
 - New megahit failure handling (#194)
 - Enforce sample wildcard constraints in Snakemake rules (#190)
 - Run megahit multithreaded (#189)

#### v2.0.2 (August 28, 2019)

 - Add implicit dependencies (samtools and bcftools) to environment file to make them explicit

#### v2.0.1 (July 24, 2019)

 - Increment Snakemake version requirement for compatibility with recent conda
 - Specify earlier megahit version to ensure compatbility with existing assembly behavior
 - Integration test improvements

#### v2.0.0 (January 22, 2019)

 - Start a project using resources directly from the SRA using `sunbeam init --data_acc [SRA ###]`. For more information, see [the docs](https://sunbeam.readthedocs.io/en/latest/usage.html#creating-a-new-project-using-data-from-sra)
 - New extension website: https://www.sunbeam-labs.org/
 - Improved documentation
 - Numerous bugfixes and optimizations

#### v1.2.1 (May 24, 2018)

 - Minor bugfixes

#### v1.2.0 (May 2, 2018)

 - Low-complexity reads are now removed by default rather than masked
 - Bug fixes related to single-end sequencing experiments
 - Documentation updates
 
#### v1.1.0 (April 8, 2018)

 - Reports include number of filtered reads per host, rather than in aggregate
 - Static binary dependency for [komplexity](https://github.com/eclarke/komplexity) for easier deployment
 - Remove max length filter for contigs
 
#### v1.0.0 (March 22, 2018)

 - First stable release!
 - Support for single-end sequencing experiments
 - Low-complexity read masking via [komplexity](https://github.com/eclarke/komplexity)
 - Support for extensions
 - Documentation on [ReadTheDocs.io](http://sunbeam.readthedocs.io)
 - Better assembler (megahit)
 - Better ORF finder (prodigal)
 - Can remove reads from any number of host/contaminant genomes
 - Semantic versioning checks
 - Integration tests and continuous deployment

-------

### Contributors

- Erik Clarke ([@eclarke](https://github.com/eclarke))
- Chunyu Zhao ([@zhaoc1](https://github.com/zhaoc1))
- Jesse Connell ([@ressy](https://github.com/ressy))
- Louis Taylor ([@louiejtaylor](https://github.com/louiejtaylor))
- Charlie Bushman ([@ulthran](https://github.com/ulthran))
- Kyle Bittinger ([@kylebittinger](https://github.com/kylebittinger))

