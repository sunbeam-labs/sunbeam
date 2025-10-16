.. _examples:

==============
Sunbeam Examples
==============

Some examples that show sunbeam in action. These are all assuming a successful sunbeam install.

Standard Run
============

You have a dataset on a big computer in your lab and you want to run sunbeam to automate qc, decontam, and assembly of metagenomic contigs. Your data is paired end and lives in a directory called ``/data``. Run:

.. code-block:: bash

    sunbeam extend https://github.com/sunbeam-labs/sbx_assembly
    sunbeam init --data_fp /data/ /projects/my_project/
    sunbeam run --profile /projects/my_project all_assembly

Once this run completes, you will have a directory called ``/projects/my_project/sunbeam_output/`` that contains all of the output from the run. Look in ``/projects/my_project/sunbeam_output/assembly/contigs/`` for the assembled contigs.

Running on a Slurm Cluster
==========================

You have a dataset on your institutions HPC (which uses slurm to manage jobs) and you want to run sunbeam to automate qc, decontam, and assembly of metagenomic contigs. Your data is paired end and lives in a directory called ``/data``. Run:

.. code-block:: bash

    sunbeam extend https://github.com/sunbeam-labs/sbx_assembly
    sunbeam init --data_fp /data/ --profile slurm /projects/my_project/
    pip install snakemake-executor-plugin-slurm

There are a number of ways to go about running this. I'll give my preferred in detail, but mention the others here:

1) Submit the main snakemake process to the cluster as well. This is generally frowned upon by HPC admins but can be a nice solution. WARNING: Recent updates to the snakemake slurm executor plugin seem to have partially/totally broken this method.
2) Use tmux/screen sessions. This is probably the best solution but I don't use it as often. Definitely recommend tmux over screen because it is more modern and handles all the environment patching better.
3) Use nohup. This is what I usually do because it is simple and works well enough.

Now you're ready to run the project, but you want to be able to logoff during the (possibly long) run, so you create a bash script called ``run_sunbeam.sh``:

.. code-block:: bash

    #!/bin/bash
    
    nohup sunbeam run --profile /projects/my_project all_assembly > log 2>&1 &

Then you submit the job:

.. code-block:: bash

    bash run_sunbeam.sh

Once this run completes, you will have a directory called ``/projects/my_project/sunbeam_output/`` that contains all of the output from the run and ``.snakemake/slurm_logs/rule_*/`` directories wherever you ran the main script from that contain logs for each job. Look in ``/projects/my_project/sunbeam_output/assembly/contigs/`` for the assembled contigs.

Using Containerized Environments
--------------------------------

Same setup, but one of the conda environments snakemake tries to make is failing to solve. You switch over to using the containerized environments and try again:

Edit /projects/my_project/config.yaml to set software-deployment-method to:

.. code-block:: yaml

    software-deployment-method: "apptainer"

And run:

.. code-block:: bash

    sunbeam run --profile /projects/my_project all_assembly

.. tip::

    Most conda environments we use are underspecified, meaning that the conda solver is left mostly to its own devices. The advantage of this is that it automates getting the most up to date versions of dependencies. The disadvantage is that sometimes the solver can't find a solution and whenever the environment does change, it risks breaking how we depend on it in sunbeam.

    Using the containerized environments guarantees that the dependencies will remain the same everytime. As long as you don't need more updated versions of the dependencies and you can run singularity or apptainer, containerization is the way to go. You can also run ``sunbeam init --data_fp ... --profile apptainer ...`` to set up the project to use containerized environments from the start.

Using the Containerized Install
===============================

You installed sunbeam via Docker, and you want to run qc and decontam then map reads back onto reference genomes. Your data is paired end and lives in a directory called ``/data`` and your reference genomes are in ``/ref_genomes``. You want the pipeline outputs in a directory called ``/projects``. To install ``sbx_mapping`` and run the pipeline, run:

.. code-block:: bash

    docker run -v /data:/data -v /ref_genomes:/ref_genomes -v /projects:/projects -it --name sunbeam sunbeamlabs/sunbeam:latest /bin/bash

    ### WITHIN THE CONTAINER ###
    sunbeam extend https://github.com/sunbeam-labs/sbx_mapping.git
    sunbeam init --data_fp /data/ /projects/my_project/
    vi /projects/my_project/config.yaml  # edit the config file to point to your reference genomes and make any other deisred changes
    sunbeam run --profile /projects/my_project all_mapping
    exit

    ### BACK OUTSIDE THE CONTAINER ###
    ls /projects/my_project/sunbeam_output/mapping/

.. tip::

    The ``-v`` flag mounts the directories from the host machine into the container. The ``-it`` flag makes the container interactive, so you can run commands within it. The ``sunbeamlabs/sunbeam:latest`` is the latest version of sunbeam, but you can also specify a version number if you want to use a specific version.

You should now be able to see all the coverage reports and other outputs from the mapping run in ``/projects/my_project/sunbeam_output/mapping/``. Note that the ``sbx_mapping`` extension was installed in this container, NOT in the image itself, so if you delete this container or start a new one, you will need to install the extension again.

.. tip::

    Including any other database in a containerized run is as simple as mounting the database directory into the container and pointing to it in the config file. For example, if you have a kraken database in ``/kraken_db`` and you want to use it in a containerized run, you would add ``-v /kraken_db:/kraken_db`` to the ``docker run`` command and then set the ``kraken_db_fp`` parameter in the config file to ``/kraken_db`` (after installing ``sbx_kraken``).

Using Slurm with the Containerized Install and Containerized Environments
-------------------------------------------------------------------------

This is a combination of all of the above. You have a dataset on your institution's HPC (which uses slurm to manage jobs) and you want to run sunbeam to automate qc, decontam, and assembly of metagenomic contigs. Your data is paired end and lives in a directory called ``/data``. You have singularity available on the cluster and you want to use containerized environments (apptainer or similar solutions also work). Run:

.. code-block:: bash

    ### FIND AND LOAD SINGULARITY MODULE IF IT ISN'T ALREADY LOADED ###
    modulefiles_list | grep singularity
    module load singularity

    export SINGULARITY_TMPDIR=/path/to/tmpdir
    export SINGULARITY_CACHEDIR=/path/to/cachedir
    singularity shell -B /data:/data,/projects:/projects docker://sunbeamlabs/sunbeam

    ### WITHIN THE CONTAINER ###
    sunbeam extend https://github.com/sunbeam-labs/sbx_assembly
    sunbeam init --data_fp /data/ --profile slurm /projects/my_project/
    pip install snakemake-executor-plugin-slurm

    ### BACK OUTSIDE THE CONTAINER ###
    vi /projects/my_project/config.yaml  # edit the config file to make any desired changes (including switching ``software-deployment-method`` to ``apptainer``)

.. note::

    Setting the ``SINGULARITY_TMPDIR`` environment variable may be necessary to avoid a bug in singularity that causes it to fail when running snakemake. The path should be a directory that is writable by the user running the container and large enough not to run out of space.

Now you're ready to run the project, but you want to be able to put the main node on the cluster as well so you can logoff, so you create a bash script called ``run_sunbeam.sh``:

.. code-block:: bash

    #!/bin/bash
    #SBATCH --time=72:00:00
    #SBATCH -n 1
    #SBATCH --mem=8G
    #SBATCH --mail-type=END,FAIL
    #SBATCH --mail-user=your_email@your_institution.edu
    #SBATCH --no-requeue
    #SBATCH --output=slurm_%x_%j.out

    # Conda env must be active for this to work
    set -x
    set -e
    singularity run -B /data:/data,/projects:/projects docker://sunbeamlabs/sunbeam sunbeam run --profile /projects/my_project all_assembly

Then you submit the job:

.. code-block:: bash

    sbatch run_sunbeam.sh

Skipping the QC and Decontamination
===================================

This time you're coming at sunbeam with a data set that you have already run QC on and removed host reads from. You want to run the assembly pipeline on this data. Your data is paired end and lives in a directory called ``/data``. Run:

.. code-block:: bash

    sunbeam extend https://github.com/sunbeam-labs/sbx_assembly
    sunbeam init --data_fp /data/ /projects/my_project/
    sunbeam run --profile /projects/my_project --skip decontam all_assembly

Once this run completes, you will have a directory called ``/projects/my_project/sunbeam_output/`` that contains all of the output from the run. Look in ``/projects/my_project/sunbeam_output/assembly/contigs/`` for the assembled contigs.

Running on AWS Batch with AWS S3 Data
======================================

COMING SOON!!!