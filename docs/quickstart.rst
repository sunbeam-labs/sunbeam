.. _quickstart:

=====================
Quickstart Guide
=====================

Installation
************

Sunbeam is a Python package that can be installed in a variety of ways.

.. tip::

   The core of Sunbeam is written in Snakemake which means it might not behave well if you're trying to integrate it into a larger pipeline. Instead consider making a Sunbeam extension for what you need.

.. tabs::

   .. tab:: pip

      From the PyPi repository (Python 3.11+ required):

      .. code-block:: shell

         python -m venv --python=python3.13 sunbeam_env/
         source sunbeam_env/bin/activate
         pip install sunbeamlib

         sunbeam -h

   .. tab:: conda

      From the BioConda channel:

      .. code-block:: shell

         conda create --name sunbeam_env -c conda-forge -c bioconda sunbeamlib

         conda activate sunbeam_env
         sunbeam -h

   .. tab:: git

      Clone and install directly from GitHub:

      .. code-block:: shell

         git clone https://github.com/sunbeam-labs/sunbeam.git
         cd sunbeam
         python -m venv --python=python3.13 sunbeam_env/
         source sunbeam_env/bin/activate
         pip install -e .[dev]

         sunbeam -h

      .. tip::

         If you're planning on doing development work on sunbeam, include the ``[dev]`` to indicate you want the extra requirements installed.
   
   .. tab:: docker

      Pull and run the Sunbeam Docker image. You'll need to have `Docker <https://docs.docker.com/get-docker/>`_ installed and running (or an alternative like Singularity or Apptainer).

      .. code-block:: shell

         docker pull sunbeamlabs/sunbeam:latest

         docker run --rm sunbeamlabs/sunbeam:latest sunbeam -h

      .. tip::

         The ``--rm`` flag removes the container it creates to run the command after it's done. This way you don't end up with a pile of dead containers on your machine. There are multiple sunbeam images available including the default which comes with prebuilt conda environments and the ``slim`` version which is smaller but requires you to build the conda environments yourself. See the `Docker Hub <https://hub.docker.com/r/sunbeamlabs/sunbeam>`_ for more information.

.. tip::

   Refer to the examples page for lots of walkthroughs of common Sunbeam use cases.

Project Initialization
**********************

Let's say your sequencing reads live in a folder called ``/sequencing/project/reads``, with one or two files per sample (for single- and paired-end sequencing, respectively). These files *must* be in gzipped FASTQ format (.fastq.gz).

Let's create a new Sunbeam project (we'll call it ``my_project``):

.. tip::

   These commands pick up where the installation instructions left off. If you're in a virtual environment, you should still be in it. If you're using Docker, you should have already run the container and be inside it.

.. tabs::

   .. tab:: Standard (Conda, local)

      Using Conda to manage worker environments and keeping all compute local:

      .. code-block:: shell

         sunbeam init my_project --data_fp /sequencing/project/reads
   
   .. tab:: Slurm

      Using Conda to manage worker environments and submitting jobs to a Slurm cluster:

      .. code-block:: shell

         pip install snakemake-executor-plugin-slurm
         sunbeam init my_project --data_fp /sequencing/project/reads --profile slurm

   .. tab:: Apptainer/Singularity

      Using Apptainer/Singularity to manage worker environments and keeping all compute local:

      .. code-block:: shell

         sunbeam init my_project --data_fp /sequencing/project/reads --profile apptainer

   .. tab:: Docker

      Using the Sunbeam Docker image to run the pipeline and keeping all compute local:

      .. code-block:: shell

         docker run --rm -v /local/path/to/data/:/data/ -v /local/path/to/outputs/:/projects/ sunbeamlabs/sunbeam:latest sunbeam init --data_fp /data/reads/ /projects/my_project

      .. tip::

         The ``-v`` flag mounts a local directory to the container. This way you can access your data and outputs from inside the container. The first ``/local/path/to/data/`` is where your data is stored on your local machine, and the second ``/local/path/to/outputs/`` is where you want the output to be saved. The ``/data/`` and ``/projects/`` are the paths inside the container that correspond to those directories.

.. tip::

   Snakemake has a number of different options for environment managers, compute services, and storage backends. See docs on executor and storage plugins for more information. And remember that you have to install the relevant plugin before you can run it.

Sunbeam will create a new folder called ``my_project`` and put three files there:

- ``config.yaml`` contains a `snakemake profile <https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles>`_ that will be used to run ``my_project``.

- ``sunbeam_config.yml`` contains all the configuration parameters for each step of the Sunbeam pipeline.

- ``samples.csv`` is a comma-separated list of samples that Sunbeam found in the given data folder, along with absolute paths to their FASTQ files.

Right now we have everything we need to do basic quality-control. However, let's go ahead and set up contaminant filtering to make things interesting.

Contaminant filtering
---------------------

Sunbeam can align your reads to an arbitrary number of contaminant sequences or host genomes and remove reads that map above a given threshold.

To use this, make a folder containing all the target sequences in FASTA format. The filenames should end in ``.fasta`` to be recognized by Sunbeam. In your ``sunbeam_config.yml`` file, edit the ``host_fp:`` line in the ``qc`` section to point to this folder.

Running the Pipeline
********************

.. tip::

   If you installed Sunbeam using Pip, you will need to have either Conda or Apptainer/Singularity installed to run the pipeline, depending on your choice of dependency manager (conda is the default).

After you've finished editing your config file, you're ready to run Sunbeam:

.. tabs::

   .. tab:: Most cases

      In most cases (Standard, Slurm, Apptainer/Singularity from the Init instructions), you can run the pipeline with:

      .. code-block:: bash

         sunbeam run --profile my_project/

   .. tab:: Docker

      If you're running Sunbeam from the Docker image, you need to be sure to mount the project directory and any database directories you want to use. Also make sure paths in your config are correct for the container, NOT your local machine.

      .. code-block:: bash

         docker run --rm -v /local/path/to/outputs/:/projects/ -v /local/path/to/blast_db/:/blast_db/ sunbeamlabs/sunbeam:latest sunbeam run --profile /projects/my_project/

      .. tip::
         
         If you're using the ``slim`` image, you will want to consider where your conda environments are stored. You could mount a local directory specifically for storing these and then point to it with ``sunbeam run --conda-prefix /conda_envs/ ...``. Or you could run ``docker run --name sunbeam ...`` without the ``--rm`` and persist the same container across runs. Or just resolve the environments on every run, which is slow and network intensive but maybe you have your reasons.

By default, this will do a lot, including trimming and quality-controlling your
reads and removing contaminant, host, and low-complexity sequences.

Viewing Results
***************

The output is stored under ``my_project/sunbeam_output``. QCed and decontaminated reads are in ``my_project/sunbeam_output/qc/decontam/``.

Extending the Pipeline
**********************

See the :ref:`extensions` page for instructions on how to add extensions to your Sunbeam project.