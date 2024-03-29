.. _quickstart:

=====================
Quickstart Guide
=====================

Installation
************

There are multiple installation methods available. We do not currently support non-Linux environments but Sunbeam has been run successfully on macOS and on Windows when using WSL.

.. tabs::

   .. tab:: tar

      On a Linux machine, download the tarball for the sunbeam version you want (``sunbeamX.X.X``) then unpack and install it.

      .. code-block:: shell

         wget https://github.com/sunbeam-labs/sunbeam/releases/latest/download/sunbeam.tar.gz
         mkdir sunbeam4.0.0
         tar -zxf sunbeam.tar.gz -C sunbeam4.0.0
         cd sunbeam4.0.0 && ./install.sh

   .. tab:: git

      On a Linux machine, download a copy of Sunbeam from our GitHub repository, and install.

      .. code-block:: shell

         git clone --branch v4.0.0 https://github.com/sunbeam-labs/sunbeam.git
         cd sunbeam
         ./install.sh

      .. tip::

         If you're planning on doing development work on sunbeam, use SSH or HTTPS with a PAT.
   
   .. tab:: docker

      On a Linux machine, you can use the Sunbeam Docker image. You'll need to have `Docker <https://docs.docker.com/get-docker/>`_ installed and running.

      .. code-block:: shell

         docker pull sunbeamlabs/sunbeam:latest
         docker run -v /local/path/to/data/:/mnt/data/ -v /local/path/to/outputs/:/mnt/projects/ -it sunbeamlabs/sunbeam:latest /bin/bash

         ### WITHIN THE CONTAINER ###
         pytest tests/  # To verify the installation
         exit

      This will drop you into a shell inside the Sunbeam Docker container. You can then run Sunbeam as you would on a normal machine.

      .. tip::

         If you're not already familiar with Docker, you may want to read up on it or use a different installation method.

This installs Sunbeam and all its dependencies, including the `Conda <https://conda.io/miniconda.html>`_ environment manager, if required. It will finish by printing instructions to continue that should look like:

.. code-block:: shell

   conda activate ENV_NAME
   pytest tests/

This runs some tests to make sure everything was installed correctly.

.. tip::

   If you've never installed Conda before, you'll need to add it to your shell's path. If you're running Bash (the most common terminal shell), the installation script should print the necessary command.

If the tests fail, check out our troubleshooting section or file an issue on our `GitHub <https://github.com/sunbeam-labs/sunbeam/issues>`_ page.

.. tip::

   Refer to the examples page for lots of walkthroughs of common Sunbeam use cases.

Setup
*****

Let's say your sequencing reads live in a folder called ``/sequencing/project/reads``, with one or two files per sample (for single- and paired-end sequencing, respectively). These files *must* be in gzipped FASTQ format.

Let's create a new Sunbeam project (we'll call it ``my_project``):

.. tabs::

   .. tab:: Standard (Conda, local)

      .. code-block:: shell

         source activate ENV_NAME
         sunbeam init my_project --data_fp /sequencing/project/reads
   
   .. tab:: Slurm

      .. code-block:: shell

         source activate ENV_NAME
         sunbeam init my_project --data_fp /sequencing/project/reads --profile slurm

   .. tab:: Apptainer/Singularity

      .. code-block:: shell

         source activate ENV_NAME
         sunbeam init my_project --data_fp /sequencing/project/reads --profile apptainer

Sunbeam will create a new folder called ``my_project`` and put three files there:

- ``config.yaml`` contains a `snakemake profile <https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles>`_ that will be used to run ``my_project``.

- ``sunbeam_config.yml`` contains all the configuration parameters for each step of the Sunbeam pipeline.

- ``samples.csv`` is a comma-separated list of samples that Sunbeam found in the given data folder, along with absolute paths to their FASTQ files.

Right now we have everything we need to do basic quality-control. However, let's go ahead and set up contaminant filtering to make things interesting.

Contaminant filtering
---------------------

Sunbeam can align your reads to an arbitrary number of contaminant sequences or host genomes and remove reads that map above a given threshold.

To use this, make a folder containing all the target sequences in FASTA format. The filenames should end in "fasta" to be recognized by Sunbeam. In your ``sunbeam_config.yml`` file, edit the ``host_fp:`` line in the ``qc`` section to point to this folder.

Running
*******

After you've finished editing your config file, you're ready to run Sunbeam:

.. code-block:: bash

   sunbeam run --profile my_project/

By default, this will do a lot, including trimming and quality-controlling your
reads and removing contaminant, host, and low-complexity sequences. Each of these steps can also be run independently by adding arguments after the ``sunbeam run`` command. See :ref:`running` for more info.

Viewing results
***************

The output is stored by default under ``my_project/sunbeam_output``. For more information on the output files and all of Sunbeam's different parts, see our full :ref:`usage`!
