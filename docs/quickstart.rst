.. _quickstart:

=====================
Quickstart Guide
=====================

.. contents::
   :depth: 2

Installation
************

On a Linux machine, download a copy of Sunbeam from our GitHub repository, and
install. We do not currently support non-Linux environments.

.. code-block:: shell

   git clone -b stable https://github.com/sunbeam-labs/sunbeam sunbeam-stable
   cd sunbeam-stable
   ./install.sh
   tests/run_tests.bash -e sunbeam

This installs Sunbeam and all its dependencies, including the `Conda
<https://conda.io/miniconda.html>`_ environment manager, if required. It then
runs some tests to make sure everything was installed correctly.

.. tip::

   If you've never installed Conda before, you'll need to add it to your shell's
   path. If you're running Bash (the most common terminal shell), the following
   command will add it to your path: ``echo 'export
   PATH=$PATH:$HOME/miniconda3/bin' > ~/.bashrc``

If you see "Tests failed", check out our :ref:`troubleshooting` section or file an issue
on our `GitHub <https://github.com/sunbeam-labs/sunbeam/issues>`_ page.

Setup
*****

Let's say your sequencing reads live in a folder called
``/sequencing/project/reads``, with one or two files per sample (for single- and
paired-end sequencing, respectively). These files *must* be in gzipped FASTQ
format.

Let's create a new Sunbeam project (we'll call it ``my_project``):

.. code-block:: shell

   source activate sunbeam
   sunbeam init my_project --data_fp /sequencing/project/reads

Sunbeam will create a new folder called ``my_project`` and put two files
there:

- ``sunbeam_config.yml`` contains all the configuration parameters for each step
  of the Sunbeam pipeline.

- ``samples.csv`` is a comma-separated list of samples that Sunbeam found the
  given data folder, along with absolute paths to their FASTQ files.

Right now we have everything we need to do basic quality-control and contig assembly. However, let's go ahead and set up contaminant filtering and some basic taxonomy databases to make things interesting.

Contaminant filtering
---------------------

Sunbeam can align your reads to an arbitrary number of contaminant sequences or
host genomes and remove reads that map above a given threshold.

To use this, make a folder containing all the target sequences in FASTA
format. The filenames should end in "fasta" to be recognized by Sunbeam. In your ``sunbeam_config.yml`` file, edit the ``host_fp:`` line in the ``qc``
section to point to this folder.

Taxonomic classification
------------------------

Sunbeam can use Kraken to assign putative taxonomic identities to your
reads. While creating a Kraken database is beyond the scope of this guide,
pre-built ones are available at the `Kraken homepage
<http://ccb.jhu.edu/software/kraken/>`_. Download or build one, then add the
path to the database under ``classify:kraken_db_fp:``.

Contig annotation
-----------------

Sunbeam can automatically BLAST your contigs against any number of
nucleotide or protein databases and summarize the top hits. Download or create
your BLAST databases, then add the paths to your config file, following the
instructions on here: :ref:`blastdbs`.

Reference mapping
-----------------

If you'd like to map the reads against a set of reference genomes of interest,
follow the same method as for the host/contaminant sequences above. Make a
folder containing FASTA files for each reference genome, then add the path to
that folder in ``mapping:genomes_fp:``.

Running
*******

After you've finished editing your config file, you're ready to run Sunbeam:

.. code-block:: bash

   sunbeam run --configfile my_project/sunbeam_config.yml

By default, this will do a lot, including trimming and quality-controlling your
reads, removing contaminant, host, and low-complexity sequences, assigning
read-level taxonomy, assembling the reads in each sample into contigs, and then
BLASTing those contigs against your databases. Each of these steps can also be run independently by adding arguments after the ``sunbeam run`` command. See :ref:`running` for more info. 

Viewing results
***************

The output is stored by default under ``my_project/sunbeam_output``. For more information on the output files and all of Sunbeam's different parts, see our full :ref:`usage`!
