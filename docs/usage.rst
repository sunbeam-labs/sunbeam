.. _usage:

==========
User Guide
==========

Requirements
============

- A relatively-recent Linux computer with more than 2Gb of RAM

We do not currently support Windows or Mac. (You can run this on
Windows using the Ubuntu [WSL](https://docs.microsoft.com/en-us/windows/wsl/about)).

.. _installation:
Installation
============

Sunbeam has two options for installation, either with git or with tar. For development work 
on sunbeam, use git. For standard usage, installing each version of sunbeam that you need 
from tarballs into separate directories is recommended (i.e. if you want versions 3 and 4 installed, 
you would repeat the tar install process below for sunbeam3.1.1 and sunbeam4.0.0 (or whatever specific 
versions you want)).

.. tabs::

   .. tab:: tar install

      On a Linux machine, download the tarball for the sunbeam version you want (``sunbeamX.X.X``) 
      then unpack and install it.

      .. code-block:: shell

         wget https://github.com/sunbeam-labs/sunbeam/releases/download/v4.0.0/sunbeam.tar.gz
         mkdir sunbeam4.0.0
         tar -zxf sunbeam4.0.0.tar.gz -C sunbeam4.0.0
         cd sunbeam4.0.0 && ./install.sh

   .. tab:: git install

      On a Linux machine, download a copy of Sunbeam from our GitHub repository, and
      install.

      .. code-block:: shell

         git clone --branch v4.0.0 https://github.com/sunbeam-labs/sunbeam.git
         cd sunbeam
         ./install.sh

      .. tip::

         If you're planning on doing development work on sunbeam, use 
         'git clone git@github.com:sunbeam-labs/sunbeam.git' instead.

The installer will check for and install the three components necessary for
Sunbeam to work. The first is `Conda <https://conda.io>`_, a system for
downloading and managing software environments. The second is the Sunbeam
environment, which will contain all the core dependencies. The third is the
Sunbeam library, which provides the necessary commands to run Sunbeam.

If you don't have Conda installed prior to this, you will need to add a line
(displayed during install) to your config file (usually in ``~/.bashrc`` or
``~/.profile``). Restart your terminal after installation for this to take
effect.

Testing
-------

We've included tests that should verify all the dependencies are
installed and Sunbeam can run properly. We strongly recommend running this after
installing or updating Sunbeam:

.. code-block:: shell

   python -m pytest tests/ -vvl

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

Sunbeam v3+ is designed to be installable separately on a system that already 
has sunbeam installed. This means multiple versions of sunbeam can be installed 
on the same machine in different repositories.

.. _uninstall:
Uninstalling or reinstalling
----------------------------

If things go awry and updating doesn't work, simply uninstall and reinstall Sunbeam.

   .. code-block:: shell

      source deactivate
      conda remove -n sunbeamX.X.X --all
      cd ../ && rm -rf sunbeam/

Then follow the installation_ instructions above.

Installing Sunbeam extensions
-----------------------------

As of version 3.0, Sunbeam extensions can be installed by running ``sunbeam extend``
followed by the URL of the extension's GitHub repo::

    sunbeam extend https://github.com/sunbeam-labs/sbx_kaiju/

For Sunbeam versions prior to 3.0, follow the legacy installation instructions on the extension to
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
  'conda env list'.

Creating a new project using local data
----------------------

We provide a utility, ``sunbeam init``, to create a new config file, profile, and sample
list for a project. The utility takes one required argument: a path to your
project folder. This folder will be created if it doesn't exist. You can also
specify the path to your gzipped fastq files, and Sunbeam will try to guess how
your samples are named, and whether they're paired.

.. code-block:: shell

   sunbeam init --data_fp /path/to/fastq/files /path/to/my_project

In this directory, a new config file and a new sample list were created (by
default named ``sunbeam_config.yml`` and ``samplelist.csv``, respectively) as well as a 
profile file (named ``config.yaml``). Edit
the config and profile files in your favorite text editor. All the keys for the config are 
described below.

.. note::

   Sunbeam will do its best to determine how your samples are named in the
   ``data_fp`` you specify. It assumes they are named something regular, like
   ``MP66_S109_L008_R1.fastq.gz`` and ``MP66_S109_L008_R2.fastq.gz``. In
   this case, the sample name would be 'MP66_S109_L008' and the read pair
   indicator would be '1' and '2'. Thus, the filename format would look like
   ``{sample}_R{rp}.fastq.gz``, where {sample} defines the sample name and
   {rp} defines the 1 or 2 in the read pair.

   If you have single-end reads, you can pass ``--single_end`` to ``sunbeam
   init`` and it will not try to identify read pairs.

   If the guessing doesn't work as expected, you can manually specify the
   filename format after the ``--format`` option in ``sunbeam init``.

   Finally, if you don't have your data ready yet, simply omit the ``--data_fp``
   option. You can create a sample list later with ``sunbeam list_samples > samples.csv``.

If some config values are always the same for all projects (e.g. paths to shared
databases), you can put these keys in a file and auto-populate your config file
with them during initialization. For instance, if you have a custom trimmomatic adapter template 
located at ``/home/user/adapter.fa``, you could have a file containing the
following called ``common_values.yml``:

.. code-block:: yaml

   qc:
     adapter_template: "/home/user/adapter.fa"

When you make a new Sunbeam project, use the ``--defaults common_values.yml`` as
part of the init command.

If you have Sunbeam extensions installed, in Sunbeam >= 3.0, the extension config
options will be automatically included in new config files generated by
``sunbeam init``.

If you want to customize options in the profile instead, you can create a custom profile 
template named ``sunbeamlib/data/custom_profile.yaml`` and fill it with whatever options you 
want included in each sunbeam run. Snakemake has a curated list of common profiles 
`here <https://github.com/Snakemake-Profiles>`_ for working with HPC platforms and job schedulers. 
A default and a slurm profile are included by default. You would use this custom profile with 
``--profile custom`` as part of the init command.

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
* ``version``: Automatically added for you by ``sunbeam init``. Ensures
  compatibility with the right version of Sunbeam.

qc
++++

* ``suffix``: the name of the subfolder to hold outputs from the
  quality-control steps
* ``leading``: (trimmomatic) remove the leading bases of a read if below this
  quality
* ``trailing``: (trimmomatic) remove the trailing bases of a read if below
  this quality
* ``slidingwindow``: (trimmomatic) the [width, avg. quality] of the sliding
  window
* ``minlength``: (trimmomatic) drop reads smaller than this length
* ``adapter_template``: (trimmomatic) path to the Illumina paired-end adaptors (templated with ``$CONDA_ENV``)
  (autofilled)
* ``fwd_adapters``: (cutadapt) custom forward adaptor sequences to remove
  using cutadapt. Replace with ``""`` to skip.
* ``rev_adapters``: (cutadapt) custom reverse adaptor sequences to remove
  using cutadapt. Replace with ``""`` to skip.
* ``cutadapt_opts``: (cutadapt) options to pass to cutadapt. Replace with ``""`` to pass no extra options.
* ``kz_threshold``: a value between 0 and 1 to determine the low-complexity boundary (1 is most stringent). Ignored if not masking low-complexity sequences.
* ``host_fp``: the path to the folder with host/contaminant genomes (ending in
  *.fasta)

classify
++++++++

  * ``suffix``: the name of the subfolder to hold outputs from the taxonomic
    classification steps

assembly
++++++++

* ``suffix``: the name of the folder to hold outputs from the assembly steps

annotation
++++++++++

* ``suffix``: the name of the folder to hold contig annotation results

.. _blastdbs:

blastdbs
++++++++

* ``root_fp``: path to a directory containing BLAST databases (if they're all in the same place)

mapping
+++++++

* ``suffix``: the name of the subfolder to create for mapping output (bam files, etc)

benchmarks
++++++++++

* ``suffix``: the name of the subfolder to create for benchmark data

logs
++++

* ``suffix``: the name of the subfolder to create for logs

.. _dbs:

Building Databases
==================

A detailed discussion on building databases for tools used by Sunbeam, while important,
is beyond the scope of this document. Please see the following resources for more details:

* `BLAST databases <https://www.ncbi.nlm.nih.gov/books/NBK279688/>`_
* `kraken databases <https://ccb.jhu.edu/software/kraken/MANUAL.html#kraken-databases>`_
* `kraken2 databases <https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual>`_

.. tip::

  These were all moved to extensions in sunbeam v4. Some vestiges remain in the main pipeline 
  for compatibility with extensions but these should be considered deprecated and will be 
  removed in future versions.

.. _running:

Running
=======

To run Sunbeam, make sure you've activated the sunbeam environment. Then run:

.. code-block:: shell

   sunbeam run --profile path/to/project/

There are many options that you can use to determine which outputs you want. By
default, if nothing is specified, this runs the entire pipeline. However, each
section is broken up into subsections that can be called individually, and will
only execute the steps necessary to get their outputs. These are specified after
the command above and consist of the following:

* ``all_qc``: basic quality control on all reads (no host read removal)
* ``all_decontam``: quality control and host read removal on all samples

To use one of these options, simply run it like so:

.. code-block:: shell

   sunbeam run --profile path/to/project/ all_qc

In addition, since Sunbeam is really just a set of `snakemake
<http://snakemake.readthedocs.io/en/latest/executable.html>`_ rules, all the
(many) snakemake options apply here as well. Some useful ones are:

* ``-n`` performs a dry run, and will just list which rules are going to be
  executed without actually doing so.
* ``-k`` allows the workflow to continue with unrelated rules if one produces an
  error (useful for malformed samples).
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
cluster nodes, you have to use the ``--cluster`` and ``--jobs`` flags. This is 
handled by using a cluster profile instead of the default. Sunbeam comes with a 
slurm profile template but you can create others or use existing ones from 
`here <https://github.com/Snakemake-Profiles>`_. Once you've initialized a 
project with a cluster profile, run it as normal:

.. code-block:: shell

   sunbeam run --profile /path/to/cluster/project/

Edit any options set in the profile as if they are snakemake command line arguments.

Outputs
=======

This section describes all the outputs from Sunbeam. Here is an example output
directory.

.. code-block:: shell

  ├ sunbeam_output
  ├ logs
	└ qc
	    ├ cleaned
	    ├ decontam
	    ├ log
	    │   ├ decontam
	    │   ├ cutadapt
	    │   └ trimmomatic
	    └ reports

Quality control
---------------

.. code-block:: shell

   	└ qc
      ├ 00_samples
      ├ 01_cutadapt
      ├ 02_trimmomatic
      ├ 03_komplexity
	    ├ cleaned
	    ├ decontam
	    ├ log
	    │   ├ decontam
	    │   ├ komplexity
	    └ reports


This   folder   contains  the   trimmed,   low-complexity   filtered  reads   in
``cleaned``. The ``decontam`` folder contains the cleaned reads that did not map
to any contaminant or host genomes. In general, most downstream steps should reference the ``decontam`` reads.

