.. _usage:

==========
User Guide
==========

.. contents::
   :depth: 3

Requirements
============

- A relatively-recent Linux computer with more than 2Gb of RAM

We do not currently support Windows or Mac. (You can run this on
Windows using the Ubuntu [WSL](https://docs.microsoft.com/en-us/windows/wsl/about)).

.. _installation:
Installation
============

Clone the stable branch of Sunbeam and run the installation script:

.. code-block:: shell

   git clone -b stable https://github.com/sunbeam-labs/sunbeam sunbeam
   cd sunbeam
   bash install.sh

The installer will check for and install the three components necessary for
Sunbeam to work. The first is `Conda <https://conda.io>`_, a system for
downloading and managing software environments. The second is the Sunbeam
environment, which will contain all the core dependencies. The third is the
Sunbeam library, which provides the necessary commands to run Sunbeam.

All of this is handled for you automatically. If Sunbeam is already installed,
you can manage versions using the manage-version.sh script, as described below 
in the updating_ section.

If you don't have Conda installed prior to this, you will need to add a line
(displayed during install) to your config file (usually in ``~/.bashrc`` or
``~/.profile``). Restart your terminal after installation for this to take
effect.

Testing
-------

We've included a test script that should verify all the dependencies are
installed and Sunbeam can run properly. We strongly recommend running this after
installing or updating Sunbeam:

.. code-block:: shell

   bash tests/run_tests.bash

If the tests fail, you should either refer to our troubleshooting_ guide or file
an issue on our `Github page <https://github.com/sunbeam-labs/sunbeam/issues>`_.

.. tip::

  You can speed up the testing process by using the environment created during 
  the install process with something like this 
  'bash tests/run_tests.bash -e SUNBEAM_ENV_NAME'. Without this argument the 
  script will create a temporary environment.

.. _updating:
Updating
--------

Sunbeam follows semantic versioning practices. In short, this means that the
version has three numbers: major, minor and patch. For instance, a version
number of 1.2.1 has 1 as the major version, 2 as the minor, and 1 as the patch.

When we update Sunbeam, if your config files and environment will work between
upgrades, we will increment the patch or minor numbers (e.g. 1.0.0 ->
1.1.0). All you need to do is the following:

.. code-block:: shell

   git pull
   ./install.sh --upgrade all

Sunbeam v3 is designed to be installable separately on a system that already 
has sunbeam 2 installed. Follow the v3 installation instructions to run both versions 
side by side.

As of v3.1.0, the manage-version.sh script can be used to install and switch between 
different versions using './manage-version.sh -s VERSION_ID'. You can see documentation on how to use this script running 
'./manage-version.sh -h'.

.. tip::

  With './manage-version.sh -s VERSION_ID', you can use a version identifier 
  (i.e. v3.1.0), dev, stable, or another branch name 
  (i.e. 342-automate-switching-between-versions-of-sunbeam). 
  './manage-version.sh -l available' will list available identifiers.

It's a good idea to re-run the tests after using this to make sure everything is working.

.. _uninstall:
Uninstalling or reinstalling
----------------------------

If things go awry and updating doesn't work, simply uninstall and reinstall Sunbeam.

   .. code-block:: shell

      source deactivate
      ./manage-version.sh -r SUNBEAM_ENV_NAME
      rm -rf sunbeam

Then follow the installation_ instructions above.

Installing Sunbeam extensions
-----------------------------

As of version 3.0, Sunbeam extensions can be installed by running ``sunbeam extend``
followed by the URL of the extension's GitHub repo::

    sunbeam extend https://github.com/sunbeam-labs/sbx_kaiju/

For Sunbeam versions prior to 3.0, follow the instructions on the extension to
install.

Setup
=====

Activating Sunbeam
------------------

Almost all commands from this point forward require us to activate the Sunbeam
conda environment:

.. code-block:: shell

   source activate SUNBEAM_ENV_NAME

You should see '(SUNBEAM_ENV_NAME)' in your prompt when you're in the environment. To leave
the environment, run ``source deactivate`` or close the terminal.

.. tip::

  You can see a list of installed sunbeam environments using the command 
  './manage-version.sh -l installed'.

Creating a new project using local data
----------------------

We provide a utility, ``sunbeam init``, to create a new config file and sample
list for a project. The utility takes one required argument: a path to your
project folder. This folder will be created if it doesn't exist. You can also
specify the path to your gzipped fastq files, and Sunbeam will try to guess how
your samples are named, and whether they're paired.

.. code-block:: shell

   sunbeam init --data_fp /path/to/fastq/files /path/to/my_project

In this directory, a new config file and a new sample list were created (by
default named ``sunbeam_config.yml`` and ``samplelist.csv``, respectively). Edit
the config file in your favorite text editor- all the keys are described below.

.. note::

   Sunbeam will do its best to determine how your samples are named in the
   ``data_fp`` you specify. It assumes they are named something regular, like
   ``MP66_S109_L008_R1_001.fastq.gz`` and ``MP66_S109_L008_R2_001.fastq.gz``. In
   this case, the sample name would be 'MP66_S109_L008' and the read pair
   indicator would be '1' and '2'. Thus, the filename format would look like
   ``{sample}_R{rp}_001.fastq.gz``, where {sample} defines the sample name and
   {rp} defines the 1 or 2 in the read pair.

   If you have single-end reads, you can pass ``--single_end`` to ``sunbeam
   init`` and it will not try to identify read pairs.

   If the guessing doesn't work as expected, you can manually specify the
   filename format after the ``--format`` option in ``sunbeam init``.

   Finally, if you don't have your data ready yet, simply omit the ``--data_fp``
   option. You can create a sample list later with ``sunbeam list_samples``.

If some config values are always the same for all projects (e.g. paths to shared
databases), you can put these keys in a file and auto-populate your config file
with them during initialization. For instance, if your Kraken databases are
located at ``/shared/kraken/standard``, you could have a file containing the
following called ``common_values.yml``:

.. code-block:: yaml

   classify:
     kraken_db_fp: "/shared/kraken/standard"

When you make a new Sunbeam project, use the ``--defaults common_values.yml`` as
part of the init command.

If you have Sunbeam extensions installed, in Sunbeam >= 3.0, the extension config
options will be automatically included in new config files generated by
``sunbeam init``.

Further usage information is available by typing ``sunbeam init --help``.

Configuration
=============

Sunbeam has lots of configuration options, but most don't need individual
attention. Below, each is described by section.

Sections
-------

all
++++

* ``root``: The root project folder, used to resolve any relative paths in the
  rest of the config file.
* ``output_fp``: Path to where the Sunbeam outputs will be stored.
* ``samplelist_fp``: Path to a comma-separated file where each row contains a
  sample name and one or two paths (if single- or paired-end) to raw gzipped
  fastq files. This can be created for you by ``sunbeam init`` or ``sunbeam
  list_samples``.
* ``paired_end``: 'true' or 'false' depending on whether you are using paired-
  or single-end reads.
* ``download_reads``: 'true' or 'false' depending on whether you are using reads
  from NCBI SRA.
* ``version``: Automatically added for you by ``sunbeam init``. Ensures
  compatibility with the right version of Sunbeam.

qc
++++

* ``suffix``: the name of the subfolder to hold outputs from the
  quality-control steps
* ``threads``: the number of threads to use for rules in this section
* ``seq_id_ending``: if your reads are named differently, a regular expression
  string defining the pattern of the suffix. For example, if your paired read
  ids are ``@D00728:28:C9W1KANXX:0/1`` and ``@D00728:28:C9W1KANXX:0/2``, this
  entry of your config file would be:
  ``seq_id_ending: "/[12]"``
* ``java_heapsize``: the memory available to Trimmomatic
* ``leading``: (trimmomatic) remove the leading bases of a read if below this
  quality
* ``trailing``: (trimmomatic) remove the trailing bases of a read if below
  this quality
* ``slidingwindow``: (trimmomatic) the [width, avg. quality] of the sliding
  window
* ``minlength``: (trimmomatic) drop reads smaller than this length
* ``adapter_template``: (trimmomatic) path to the Illumina paired-end adaptors (templated with ``$CONDA_ENV``)
  (autofilled)
* ``fwd_adaptors``: (cutadapt) custom forward adaptor sequences to remove
  using cutadapt. Replace with ``""`` to skip.
* ``rev_adaptors``: (cutadapt) custom reverse adaptor sequences to remove
  using cutadapt. Replace with ``""`` to skip.
* ``mask_low_complexity``: [true/false] mask low-complexity sequences with Ns
* ``kz_threshold``: a value between 0 and 1 to determine the low-complexity boundary (1 is most stringent). Ignored if not masking low-complexity sequences.
* ``kz_window``: window size to use (in bp) for local complexity
  assessment. Ignored if not masking low-complexity sequences.
* ``pct_id``: (decontaminate) minimum percent identity to host genome to
  consider match
* ``frac``: (decontaminate) minimum fraction of the read that must align to
  consider match
* ``host_fp``: the path to the folder with host/contaminant genomes (ending in
  *.fasta)


classify
++++++++

  * ``suffix``: the name of the subfolder to hold outputs from the taxonomic
    classification steps
  * ``threads``: threads to use for Kraken
  * ``kraken_db_fp``: path to Kraken database


assembly
++++++++

* ``suffix``: the name of the folder to hold outputs from the assembly steps
* ``min_len``: the minimum contig length to keep
* ``threads``: threads to use for the MEGAHIT assembler

annotation
++++++++++

* ``suffix``: the name of the folder to hold contig annotation results
* ``min_contig_length``: minimum length of contig to annotate (shorter contigs are skipped)
* ``circular_kmin``: smallest length of kmers used to search for circularity
* ``circular_kmax``: longest length of kmers used to search for circularity
* ``circular_min_length``: smallest length of contig to check for circularity

blast
+++++

* ``threads``: number of threads provided to all BLAST programs

.. _blastdbs:

blastdbs
++++++++

* ``root_fp``: path to a directory containing BLAST databases (if they're all in the same place)
* ``nucleotide``: the section to define any nucleotide BLAST databases (see tip below for syntax)
* ``protein``: the section to define any protein BLAST databases (see tip below)

  .. tip::

     The structure for this section allows you to specify arbitrary numbers of
     BLAST databases of either type. For example, if you had a local copy of nt
     and a couple of custom protein databases, your section here would look like
     this (assuming they're all in the same parent directory):

     .. code-block:: yaml

	blastdbs:
          root_fp: "/local/blast_databases"
	  nucleotide:
	    nt: "nt/nt"
	  protein:
	    vfdb: "virulence_factors/virdb"
	    card: "/some/other/path/card_db/card"

     This tells Sunbeam you have three BLAST databases, two of which live in
     ``/local/blast_databases`` and a third that lives in
     ``/some/other/path``. It will run nucleotide blast on the nucleotide
     databases and BLASTX and BLASTP on the protein databases.

mapping
+++++++

* ``suffix``: the name of the subfolder to create for mapping output (bam files, etc)
* ``genomes_fp``: path to a directory with an arbitrary number of target genomes
  upon which to map reads. Genomes should be in FASTA format, and Sunbeam will
  create the indexes if necessary.
* ``threads``: number of threads to use for alignment to the target genomes
* ``samtools_opts``: a string added to the ``samtools view`` command during
  mapping. This is a good place to add '-F4' to keep only mapped reads and
  decrease the space these files occupy.

download
++++++++
* ``suffix``: the name of the subfolder to create for download output (fastq.gz files)
* ``threads``: number of threads to use for downloading (too many at once may make NCBI unhappy)

.. _dbs:

Building Databases
==================

A detailed discussion on building databases for tools used by Sunbeam, while important,
is beyond the scope of this document. Please see the following resources for more details:

* `BLAST databases <https://www.ncbi.nlm.nih.gov/books/NBK279688/>`_
* `kraken databases <https://ccb.jhu.edu/software/kraken/MANUAL.html#kraken-databases>`_
* `kraken2 databases <https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual>`_ (used in Sunbeam v3.0 and higher)

.. _running:

Running
=======

To run Sunbeam, make sure you've activated the sunbeam environment. Then run:

.. code-block:: shell

   sunbeam run --configfile ~/path/to/config.yml

There are many options that you can use to determine which outputs you want. By
default, if nothing is specified, this runs the entire pipeline. However, each
section is broken up into subsections that can be called individually, and will
only execute the steps necessary to get their outputs. These are specified after
the command above and consist of the following:

* ``all_qc``: basic quality control on all reads (no host read removal)
* ``all_decontam``: quality control and host read removal on all samples
* ``all_mapping``: align reads to target genomes
* ``all_classify``: classify taxonomic provenance of all qc'd, decontaminated
  reads
* ``all_assembly``: build contigs from all qc'd, decontaminated reads
* ``all_annotate``: annotate contigs using defined BLAST databases

To use one of these options, simply run it like so:

.. code-block:: shell

   sunbeam run -- --configfile ~/path/to/config.yml all_classify

In addition, since Sunbeam is really just a set of `snakemake
<http://snakemake.readthedocs.io/en/latest/executable.html>`_ rules, all the
(many) snakemake options apply here as well. Some useful ones are:

* ``-n`` performs a dry run, and will just list which rules are going to be
  executed without actually doing so.
* ``-k`` allows the workflow to continue with unrelated rules if one produces an
  error (useful for malformed samples, which can also be added to the
  ``exclude`` config option).
* ``-p`` prints the actual shell command executed for each rule, which is very
  helpful for debugging purposes.
* ``--cores`` specifies the total number of cores used by Sunbeam. For example,
  if you run Sunbeam with ``--cores 100`` and each rule/processing step uses
  20 threads, it will run 5 rules at once.

.. _cluster:

Cluster options
---------------

Sunbeam inherits its cluster abilities from Snakemake. There's nothing special
about installing Sunbeam on a cluster, but in order to distribute work to
cluster nodes, you have to use the ``--cluster`` and ``--jobs`` flags. For
example, if we wanted each rule to run on a 12-thread node, and a max of 100
rules executing in parallel, we would use the following command on our cluster:

.. code-block:: shell

   sunbeam run -- --configfile ~/path/to/config.yml --cluster "bsub -n 12" -j 100 -w 90

The ``-w 90`` flag is provided to account for filesystem latency that often
causes issues on clusters. It asks Snakemake to wait for 90 seconds before
complaining that an expected output file is missing.

Outputs
=======

This section describes all the outputs from Sunbeam. Here is an example output
directory, where we had two samples (sample1 and sample2), two BLAST
databases, one nucleotide ('bacteria') and one protein ('card').

.. code-block:: shell

   sunbeam_output
	├ annotation
	│   ├ blastn
	│   │   └ bacteria
	│   │       └ contig
	│   ├ blastp
	│   │   └ card
	│   │       └ prodigal
	│   ├ blastx
	│   │   └ card
	│   │       └ prodigal
	│   ├ genes
	│   │   └ prodigal
	│   │       └ log
	│   └ summary
	├ assembly
	│   ├ contigs
	├ classify
	│   └ kraken
	│       └ raw
	├ mapping
   	│   └ genome1
	└ qc
	    ├ cleaned
	    ├ decontam
	    ├ log
	    │   ├ decontam
	    │   ├ cutadapt
	    │   └ trimmomatic
	    └ reports

In order of appearance, the folders contain the following:

Contig annotation
-----------------

.. code-block:: shell

   sunbeam_output
	├ annotation
	│   ├ blastn
	│   │   └ bacteria
	│   │       └ contig
	│   ├ blastp
	│   │   └ card
	│   │       └ prodigal
	│   ├ blastx
	│   │   └ card
	│   │       └ prodigal
	│   ├ genes
	│   │   └ prodigal
	│   │       └ log
	│   └ summary

This contains the BLAST/Diamond results in blast tabular format from the assembled contigs. ``blastn``
contains the results from directly BLASTing the contig nucleotide sequences
against the nucleotide databases. ``blastp`` and ``blastx`` use genes identified
by the ORF finding program Prodigal to search for hits in the protein databases.

The genes found from Prodigal are available in the ``genes`` folder.

Finally, the ``summary`` folder contains an aggregated report of the number and
types of hits of each contig against the BLAST databases, as well as length and
circularity.

Contig assembly
---------------

.. code-block:: shell

   	├ assembly
	│   ├ contigs


This contains the assembled contigs for each sample under 'contigs'.

Taxonomic classification
------------------------

.. code-block:: shell

	├ classify
	│   └ kraken
	│       └ raw

This contains the taxonomic outputs from Kraken, both the raw output as well as
summarized results. The primary output file is ``all_samples.tsv``, which is a
BIOM-style format with samples as columns and taxonomy IDs as rows, and number
of reads assigned to each in each cell.

Alignment to genomes
--------------------

.. code-block:: shell

   	├ mapping
   	│   └ genome1


Alignment files (in BAM format) to each target genome are contained in
subfolders named for the genome, such as 'genome1'.

Quality control
---------------

.. code-block:: shell

   	└ qc
	    ├ cleaned
	    ├ decontam
	    ├ log
	    │   ├ decontam
	    │   ├ cutadapt
	    │   └ trimmomatic
	    └ reports


This   folder   contains  the   trimmed,   low-complexity   filtered  reads   in
``cleaned``. The ``decontam`` folder contains the cleaned reads that did not map
to any contaminant or host genomes. In general, most downstream steps should reference the ``decontam`` reads.

