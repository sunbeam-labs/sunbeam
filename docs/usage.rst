.. _usage:

==========
User Guide
==========

.. contents::
   :depth: 3

.. _installation:
Installation
============

Clone the stable branch of Sunbeam and run the installation script:

.. code-block:: shell

   git clone -b stable https://github.com/eclarke/sunbeam sunbeam-stable
   cd sunbeam-stable
   bash install.sh

Testing
-------

We've included a test script that should verify all the dependencies are
installed and Sunbeam can run properly. We strongly recommend running this after
installing or updating Sunbeam:

.. code-block:: shell

   bash tests/test.sh

If the test exits with anything besides 0, something isn't working and you
should either refer to our troubleshooting_ guide or file an issue on our
`Github page <https://github.com/eclarke/sunbeam/issues>`_.

.. _updating:
Updating
--------

Sunbeam follows semantic versioning. In short, this means that the version has
three numbers: x.y.z (e.g. 1.0.2). In between major updates (where x changes),
the config file and dependencies remain the same and Sunbeam can be updated in
place. When x changes, we suggest a complete reinstallation because old config
files may not work anymore.

1. Minor version upgrades

   To update between major versions, simply pull the most recent changes into
   the ``sunbeam-stable`` directory using ``git pull``. We suggest running the
   tests after this step to verify the update.

.. _uninstall:
2. Major version upgrades
   
   To upgrade versions, we suggest uninstalling and reinstalling Sunbeam. To
   uninstall, run:

   .. code-block:: shell

      conda env remove sunbeam
      rm -rf sunbeam-stable

   Then follow the installation_ instructions above.

Setup
=====

Activating Sunbeam
------------------

Almost all commands from this point forward require us to activate the Sunbeam
conda environment:

.. code-block:: shell

   source activate sunbeam

You should see '(sunbeam)' in your prompt when you're in the environment. To leave
the environment, run ``source deactivate`` or close the terminal.

Creating a new project
----------------------

We provide a utility, ``sunbeam_init``, to create a new config file for a
project. The utility takes one required argument: a path to your project
folder. This folder may be empty, or contain a subfolder with your sequencing
data. 

.. code-block:: shell

   mkdir ~/my_project
   sunbeam_init ~/my_project > ~/my_project/sunbeam_config.yml
   
We now have a config file in that directory. If you're a member of the Bushman Lab or PennCHOP group, there are defaults available for you depending on what server you're running on. To use these, pass the ``--server`` option along with the server name. For instance, if I'm running on microb120:

.. code-block:: shell

   mkdir ~/my_project
   sunbeam_init --server microb120 ~/my_project > ~/my_project/sunbeam_config.yml


Configuration
=============

Sunbeam has lots of configuration options, but most don't need individual attention. Below, each is described by section.

Sections
-------

all
++++

* ``root``: The root project folder, used to resolve any relative paths in the
  rest of the config file.
* ``data_fp``: The path to the raw, gzipped fastq sequence files.
* ``filename_fmt``: This defines how to find the sample and read pairing
  in your samples' filenames.

  .. tip::
     If your files are in pairs like ``MP66_S109_L008_R1_001.fastq.gz``
     and ``MP66_S109_L008_R2_001.fastq.gz``, the sample name would be
     'MP66_S109_L008' and the read pair (rp) would be 'R1' and 'R2'. Thus, the
     ``filename_fmt`` would be ``{sample}_{rp}_001.fastq.gz``.

* ``samplelist_fp``: The path to a file with list of sample names (one per
  line) to work on instead of finding them in the ``data_fp`` directory. This
  is useful for only working on certain samples in a folder.
* ``subcores``: currently ignored.
* ``exclude``: A list, specified using sample names in quotes between the
  square brackets, of samples to ignore. This is useful when a sample is
  causing an error downstream and you want to skip it. For example:
  
  .. code-block:: yaml
		    
     exclude: ['bad_sample1', 'bad_sample2']

qc
++++

* ``suffix``: the name of the subfolder to hold outputs from the
  quality-control steps
* ``threads``: the number of threads to use for rules in this section
* ``java_heapsize``: the memory available to Trimmomatic
* ``leading``: (trimmomatic) remove the leading bases of a read if below this
  quality
* ``trailing``: (trimmomatic) remove the trailing bases of a read if below
  this quality
* ``slidingwindow``: (trimmomatic) the [width, avg. quality] of the sliding
  window
* ``minlength``: (trimmomatic) drop reads smaller than this length
* ``adapter_fp``: (trimmomatic) path to the Illumina paired-end adaptors
  (autofilled)
* ``fwd_adaptors``: (cutadapt) custom forward adaptor sequences to remove
  using cutadapt. Replace with "" to skip.
* ``rev_adaptors``: (cutadapt) custom reverse adaptor sequences to remove
  using cutadapt. Replace with "" to skip.
* ``pct_id``: (decontaminate) minimum percent identity to host genome to
  consider match
* ``frac``: (decontaminate) minimum fraction of the read that must align to
  consider match
* ``keep_sam``: (decontaminate) keep SAM file of host read alignment for
  debuggging
* ``method``: (decontaminate) use either BWA or BowTie for alignment
* ``human_genome_fp``: The path to the host genome for host read
  removal. Despite the name, this doesn't have to be a human genome.
* ``phix_genome_fp``: The path to the PhiX genome for PhiX removal.

classify
++++++++

  * ``suffix``: the name of the subfolder to hold outputs from the taxonomic
    classification steps
  * ``threads``: threads to use for Kraken
  * ``kraken_db_fp``: path to Kraken database
  * ``taxa_db_fp``: currently ignored

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
* ``keep_unaligned``: whether or not to keep unaligned reads


Running
=======

To run Sunbeam, make sure you've activated the sunbeam environment and are in the sunbeam folder. Then run:

.. code-block:: shell

   snakemake --configfile ~/path/to/config.yml

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

   snakemake --configfile ~/path/to/config.yml all_classify

In addition, since Sunbeam is really just a set of `snakemake <http://snakemake.readthedocs.io/en/latest/executable.html>`_ rules, all the
(many) snakemake options apply here as well. Some useful ones are:

* ``-n`` performs a dry run, and will just list which rules are going to be
  executed without actually doing so.
* ``-k`` allows the workflow to continue with unrelated rules if one produces an
  error (useful for malformed samples, which can also be added to the
  ``exclude`` config option).
* ``-p`` prints the actual shell command executed for each rule, which is very
  helpful for debugging purposes.

.. _cluster:
Cluster options
---------------

Sunbeam inherits its cluster abilities from Snakemake. There's nothing special
about installing Sunbeam on a cluster, but in order to distribute work to
cluster nodes, you have to use the ``--cluster`` and ``--jobs`` flags. For
example, if we wanted each rule to run on a 12-thread node, and a max of 100
rules executing in parallel, we would use the following command on our cluster:

.. code-block:: shell

   snakemake --configfile ~/path/to/config.yml --cluster "bsub -n 12" -j 100 -w 90

The ``-w 90`` flag is provided to account for filesystem latency that often
causes issues on clusters. It asks Snakemake to wait for 90 seconds before
complaining that an expected output file is missing.


Outputs
=======

This section describes all the outputs from Sunbeam. Here is an example output
directory, where we had two samples (sample1 and sample2), and two BLAST
databases, one nucleotide ('bacteria') and one protein ('card').

.. code-block:: shell

   sunbeam_output
	├── annotation
	│   ├── blastn
	│   │   └── bacteria
	│   │       └── contig
	│   ├── blastp
	│   │   └── card
	│   │       └── prodigal
	│   ├── blastx
	│   │   └── card
	│   │       └── prodigal
	│   ├── genes
	│   │   └── prodigal
	│   │       └── log
	│   └── summary
	├── assembly
	│   ├── sample1_assembly
	│   ├── sample2_assembly
	│   ├── log
	│   │   ├── cap3
	│   │   └── vsearch
	├── classify
	│   └── kraken
	│       └── raw
	├── mapping
	└── qc
	    ├── cutadapt
	    ├── decontam
	    ├── decontam-human
	    ├── decontam-phix
	    ├── log
	    │   ├── decontam
	    │   ├── decontam-human
	    │   └── trimmomatic
	    ├── paired
	    └── unpaired

In order of appearance, the folders contain the following:

Contig annotation
-----------------

.. code-block:: shell

   sunbeam_output
	├── annotation
	│   ├── blastn
	│   │   └── bacteria
	│   │       └── contig
	│   ├── blastp
	│   │   └── card
	│   │       └── prodigal
	│   ├── blastx
	│   │   └── card
	│   │       └── prodigal
	│   ├── genes
	│   │   └── prodigal
	│   │       └── log
	│   └── summary
   
This contains the BLAST results in XML from the assembled contigs. ``blastn``
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

	├── assembly
	│   ├── sample1_assembly
	│   ├── sample2_assembly
	│   ├── log
	│   │   ├── cap3
	│   │   └── vsearch

This contains the assembled contigs for each sample in its own folder under [samplename]_assembly.

Taxonomic classification
------------------------

.. code-block:: shell
   
	├── classify
	│   └── kraken
	│       └── raw

This contains the taxonomic outputs from Kraken, both the raw output as well as
summarized results. The primary output file is ``all_samples.tsv``, which is a
BIOM-style format with samples as columns and taxonomy IDs as rows, and number
of reads assigned to each in each cell.

Alignment to genomes
--------------------

.. code-block:: shell
   
	├── mapping

Right now this contains all the output files (.bam) for the mapping of reads
back to target genomes. We plan on breaking down the output into subfolders for
a more organized structure soon.

Quality control
---------------

.. code-block:: shell
   
	└── qc
	    ├── cutadapt
	    ├── decontam
	    ├── decontam-human
	    ├── decontam-phix
	    ├── log
	    │   ├── decontam
	    │   ├── decontam-human
	    │   └── trimmomatic
	    ├── paired
	    └── unpaired


This folder contains paired, non-host-removed reads in ``paired`` (and unpaired
in ``unpaired``). ``decontam`` contains the final output of both the host and
phix removal steps.
	

.. _troubleshooting:
Troubleshooting
===============
