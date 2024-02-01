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

    # Conda env must be activate for this to work
    set -x
    set -e
    sunbeam run --profile /projects/my_project all_assembly

Then you submit the job:

.. code-block:: bash

    sbatch run_sunbeam.sh

Once this run completes, you will have a directory called ``/projects/my_project/sunbeam_output/`` that contains all of the output from the run and ``slurm_*`` files wherever you ran the main script from that contains logs for each job. Look in ``/projects/my_project/sunbeam_output/assembly/contigs/`` for the assembled contigs.

Using Containerized Environments
--------------------------------

Same setup, but one of the conda environments snakemake tries to make is failing to solve. You switch over to using the containerized environments and try again:

Edit /projects/my_project/config.yaml to set software-deployment-method to:

.. code-block:: yaml

    software-deployment-method:
        - "conda"
        - "apptainer"

And run:

.. code-block:: bash

    sunbeam run --profile /projects/my_project all_assembly

.. tip::

    Most conda environments we use are underspecified, meaning that the conda solver is left mostly to its own devices. The advantage of this is that it automates getting the most up to date versions of dependencies. The disadvantage is that sometimes the solver can't find a solution and whenever the environment does change, it risks breaking how we depend on it in sunbeam.

    Using the containerized environments guarantees that the dependencies will remain the same everytime. As long as you don't need more updated versions of the dependencies and you can run singularity or apptainer, containerization is the way to go.