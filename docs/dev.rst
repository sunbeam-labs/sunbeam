.. _dev:

===
Dev
===

Getting involved with developing Sunbeam can be a little daunting at first. This doc will try to break down the constituent parts from a developer's perspective. For starters, check out the structure_ doc to get a sense of how the code is organized.

sunbeam (Python lib)
====================

Sunbeam is a Python package designed to help facilitate reproducible bioinformatics workflows. The core of the package is located in the ``sunbeam/`` directory. This is where the main code for Sunbeam lives, and it is organized into subdirectories based on functionality.

.. autofunction:: sunbeam.get_docker_str
.. autofunction:: sunbeam.get_ext_path

bfx
---

The ``bfx/`` directory contains bioinformatics utilties, used in the pipeline for transformation and reporting.

.. automodule:: sunbeam.bfx
    :members:

    :autofunction:: sunbeam.bfx.decontam.get_mapped_reads

configs
-------

The ``configs/`` directory contains sample configuration files for the pipeline. These are used by ``sunbeam init`` to set up a project and define the parameters for the workflow.

extensions
----------

This is the default location for extensions, although it can be configured by setting ``$SUNBEAM_EXTENSIONS``.

project
-------

The ``project/`` directory contains project management utilities, used for initializing and managing Sunbeam projects.

.. automodule:: sunbeam.project
    :members:

scripts
-------

The ``scripts/`` directory contains scripts for running the pipeline and managing the workflow.

workflow
--------

The core of the work done by Sunbeam is handled by the Snakemake workflow, located in ``workflow/``. Once a project is setup properly, the workflow can be run with all the benefits of Snakemake. The core of the workflow is defined in ``workflow/Snakefile``. Reference the Snakemake docs for help understanding Snakemake things better; they're very good. From this core Snakefile, we import more ``.smk`` files from ``workflow/rules/`` and ``extensions/sbx_*/``.

Important Variables (Python)
----------------------------

Variables defined in the Python library and can optionally be imported into Snakemake:

- ``__version__``: str - The version of Sunbeam being run.
- ``EXTENSIONS_DIR``: () -> Path - A function that lazily loads the Path to the extensions directory.
- ``WORKFLOW_DIR``: Path - The Path to the workflow directory.
- ``CONFIGS_DIR``: Path - The Path to the configs directory.

Important Variables (Snakemake)
-------------------------------

Variables defined in the main Snakefile can be accessed throughout the workflow. Some important variables include:

- ``Samples``: Dict[str, Dict[str, str]] - A dictionary where keys are sample names and values are dictionaries of read pairs mapping to file paths (``Samples[sample] = {"1": r1, "2": r2}``).
- ``Pairs``: List[str] - Either ``["1", "2"]`` or ``["1"]`` depending on if the project is paired end or not.
- ``Cfg``: Dict[str, Dict[str, str]] - The YAML config converted into dictionary form.
- ``MIN_MEM_MB``: int - A minimum value of the number of megabytes of memory to request for each job. This will only apply for jobs that rely on Sunbeam to guess their memory requirements.
- ``MIN_RUNTIME``: int - A minimum value of the number of minutes to request for each job. This will only apply for jobs that rely on Sunbeam to guess their runtime requirements.
- ``HostGenomes``: List[str] - A list of host genomes that are used for decontaminating reads.
- ``HostGenomeFiles``: List[str] - A list of files with host genomes that are used for decontaminating reads (not to be confused with ``sbx_mapping``'s ``GenomeFiles`` variable, which it uses to track reference genome files).
- ``QC_FP``: Path - The Path to the project's quality control output directory.
- ``BENCHMARK_FP``: Path - The Path to the project's benchmarking output directory.
- ``LOG_FP``: Path - The Path to the project's log output directory.

Environment Variables
---------------------

- ``SUNBEAM_EXTS_INCLUDE``: str - If set, will include the given extension in the workflow (and exclude the rest). This is useful for testing individual extensions. Can also be set in the args to ``sunbeam run``.
- ``SUNBEAM_EXTS_EXCLUDE``: str - If set, will exclude the given extension from the workflow. This is useful for when namespaces between extensions collide (same rule name multiple times). Can also be set in the args to ``sunbeam run``.
- ``SUNBEAM_SKIP``: str - If set, will skip either 'qc' or 'decontam'. Can also be set in the args to ``sunbeam run``.
- ``SUNBEAM_DOCKER_TAG``: str - If set, will use the given tag for the Docker image instead of the default. Can also be set in the args to ``sunbeam run``.
- ``SUNBEAM_MIN_MEM_MB``: int - If set, will override the default minimum memory value.
- ``SUNBEAM_MIN_RUNTIME``: int - If set, will override the default minimum runtime value.
- ``SUNBEAM_NO_ADAPTER``: bool - If set, will not check that the adapter template file exists.

tests
=====

All tests are located in the ``tests/`` directory. The tests are run with pytest, and the tests are organized into subdirectories based on the module they are testing.

.github
=======

The ``.github/`` directory contains the configuration for GitHub Actions, which are used to run the tests on every push to the repository and manage releases. The configuration is in ``.github/workflows/``.

docs
====

The ``docs/`` directory contains the documentation for Sunbeam. The documentation is written in reStructuredText and is built with Sphinx.
