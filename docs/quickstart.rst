.. _quickstart:

=====================
Quickstart Guide
=====================

.. contents::
   :depth: 2
.. tabs::

   .. code-tab:: c

      Apples are green, or sometimes red.

   .. code-tab:: py

      Pears are green.

   .. code-tab:: java

      Oranges are orange.
Installation
************

There are two installation methods available, installing via git or via tar. We do not currently support non-Linux environments.

.. tabs::

   .. tab:: tar install

      On a Linux machine, download the tarball for the sunbeam version you want (``sunbeamX.X.X``) 
      then unpack and install it.

      .. code-block:: shell

         wget https://github.com/sunbeam-labs/sunbeam/archive/refs/tags/sunbeam4.0.0.tar.gz
         mkdir sunbeam4.0.0
         tar -zxf sunbeam4.0.0.tar.gz -C sunbeam4.0.0
         cd sunbeam4.0.0 && ./install.sh

   .. tab:: git install

      On a Linux machine, download a copy of Sunbeam from our GitHub repository, and
      install.

      .. code-block:: shell

         git clone -b stable https://github.com/sunbeam-labs/sunbeam.git
         cd sunbeam
         ./install.sh

      .. tip::

         If you're planning on doing development work on sunbeam, use 
         'git clone -b stable git@github.com:sunbeam-labs/sunbeam.git' instead.

This installs Sunbeam and all its dependencies, including the `Conda
<https://conda.io/miniconda.html>`_ environment manager, if required. It will finish 
by printing instructions to continue that should look like:

.. code-block:: shell

   conda activate ENV_NAME
   tests/run_tests.bash -e ENV_NAME

This runs some tests to make sure everything was installed correctly.

.. tip::

   If you've never installed Conda before, you'll need to add it to your shell's
   path. If you're running Bash (the most common terminal shell), the installation 
   script should print the necessary command.

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

   source activate ENV_NAME
   sunbeam init my_project --data_fp /sequencing/project/reads

Sunbeam will create a new folder called ``my_project`` and put three files
there:

- ``config.yaml`` contains a `snakemake profile<https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles>`_ that will be used to run ``my_project``.

- ``sunbeam_config.yml`` contains all the configuration parameters for each step
  of the Sunbeam pipeline.

- ``samples.csv`` is a comma-separated list of samples that Sunbeam found the
  given data folder, along with absolute paths to their FASTQ files.

Right now we have everything we need to do basic quality-control and contig assembly. However, let's go ahead and set up contaminant filtering to make things interesting.

Contaminant filtering
---------------------

Sunbeam can align your reads to an arbitrary number of contaminant sequences or
host genomes and remove reads that map above a given threshold.

To use this, make a folder containing all the target sequences in FASTA
format. The filenames should end in "fasta" to be recognized by Sunbeam. In your ``sunbeam_config.yml`` file, edit the ``host_fp:`` line in the ``qc``
section to point to this folder.

Running
*******

After you've finished editing your config file, you're ready to run Sunbeam:

.. code-block:: bash

   sunbeam run --profile my_project/

By default, this will do a lot, including trimming and quality-controlling your
reads, removing contaminant, host, and low-complexity sequences, and assembling the reads in each sample into contigs. Each of these steps can also be run independently by adding arguments after the ``sunbeam run`` command. See :ref:`running` for more info. 

Viewing results
***************

The output is stored by default under ``my_project/sunbeam_output``. For more information on the output files and all of Sunbeam's different parts, see our full :ref:`usage`!
