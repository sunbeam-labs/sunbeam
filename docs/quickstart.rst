.. _quickstart:

===============
Getting Started
===============

.. contents::
   :depth: 2

Installation
************

The easiest way to get an up-to-date version of Sunbeam is to clone the latest
stable release from our source repository (other installation methods are
described in our :ref:`install` guide):

.. code-block:: shell

   git clone -b stable https://github.com/eclarke/sunbeam
   cd sunbeam
   bash install.sh

This downloads Sunbeam and installs all the dependencies. If you've never
installed Conda before, you will need to add a link to it to your profile to
make sure we can use it later. The install script should output instructions on
how to do that. Next, run:

.. code-block:: shell

   bash tests/test.sh && echo "Tests succeeded" || echo "Tests failed"

This runs the test script to make sure everything works. If you see "Tests
failed" see the :ref:`install` guide for troubleshooting or file an issue on our
`GitHub <https://github.com/eclarke/sunbeam/issues>`_ page.


Setup
*****

We need to set a few things up before we can start. Sunbeam takes a directory
full of paired, gzipped fastq files that you get from Illumina's
sequencer. Let's assume those files live in ``~/my_project/data_files``.

To set up a new Sunbeam config file, activate the Sunbeam environment and run
``sunbeam init``:

.. code-block:: shell
   
   source activate sunbeam
   sunbeam init ~/my_project > ~/my_project/config.yml

This writes a new ``config.yml`` file into your project directory, pre-filled
with some convenient defaults. However, there are some things we'll need to edit
by hand.

Sample filenames
----------------
The first thing is describing to Sunbeam what your files look like. If
we take a peek inside our ``data_files`` directory, we might see something like
this:

.. code-block:: shell

   $ ls ~/my_project/data_files
   MP66_S109_L008_R1_001.fastq.gz   VS_B3_S173_L004_R2_001.fastq.gz
   MP66_S109_L008_R2_001.fastq.gz   VS_B3_S173_L005_R1_001.fastq.gz
   MP67_S110_L001_R1_001.fastq.gz   VS_B3_S173_L005_R2_001.fastq.gz
   MP67_S110_L001_R2_001.fastq.gz   VS_B3_S173_L006_R1_001.fastq.gz
   MP67_S110_L002_R1_001.fastq.gz   VS_B3_S173_L006_R2_001.fastq.gz
   ...

In this example, our samples are the first part of this filename
(e.g. ``MP66_S109_L008``) and the read pair is the second (``R1`` or
``R2``). The rest is the same between all the files, so if we were to describe a
pattern it would look like ``{sample}_{rp}_001.fastq.gz``.

We will input this filename format into the config file now. Open the config
file in your favorite editor and change the value after ``filename_fmt:`` to
whatever best describes your samples.

Data directory
--------------

The value of ``data_fp`` points to your raw, gzipped fastq files that we described above. Edit this if the directory isn't named ``data_files``. If you use a relative path here, it will be referenced relative to the ``root`` directory (this is true for all paths in the config file).

Host genomes
------------

Sunbeam decontaminates host and PhiX reads by default. Input paths to the host and phix genome under the ``qc: human_genome_fp:`` and ``qc: phix_genome_fp:``, respectively.

.. note:: The host genome doesn't have to be human, despite the name. This should be fixed in an upcoming release.

Databases
---------

Sunbeam can use Kraken and BLAST databases to perform taxonomic and functional characterization of the reads and/or contigs after quality control and contig assembly.

To add a Kraken database, provide a path to it at ``classify: kraken_fp:``.

To add one or more BLAST databases, provide paths to them under ``blastdbs:``. You have to specify if the databases are protein or nucleotide databases; for instance, if I were using the CARD protein database and the NT nucleotide database, that section of the config file would look like this:

.. code-block:: yaml

   blastdbs:
     root_fp: "/local/blast_databases"
     nucleotide:
	nt: "nt/nt"
     protein:
	card: "card_db/card"

This tells Sunbeam that our databases are in ``/local/blast_databases/nt/nt`` and ``/local/blast_databases/card_db/card``, and that they are respectively named 'nt' and 'card'. 


Running
*******

After you've saved your edits to the config file, try a dry run of Sunbeam to see if everything looks okay. Make sure you're in the ``sunbeam`` folder and have activated Sunbeam using ``source activate sunbeam``, then run:

.. code-block:: bash

   snakemake --configfile ~/my_project/config.yml -n

The ``-n`` specifies a dry run, so Sunbeam will verify all the paths in the config file to verify they exist and complain about any errors that may occur.

If everything looks good, you can start Sunbeam by removing the ``-n``
option. At this point, the operation of Sunbeam is done using Snakemake, and all
of the (many) options available to Snakemake apply here.

Viewing results
***************

The results of Sunbeam are in the folder dictated by ``output_fp`` in the config file. The resulting folder structure contains various outputs in subdirectories like 'qc' (host-filtered and adaptor-trimmed sequences), 'classify' (Kraken classification results) and 'assembly' (reads assembled into contigs).

