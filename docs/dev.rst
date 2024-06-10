.. _dev:

====
Dev
====

Getting involved with developing Sunbeam can be a little daunting at first. This doc will try to break down the constituent parts from a developer's perspective. For starters, check out the structure_ doc to get a sense of how the code is organized.

sunbeamlib
==========

The core of Sunbeam's configuration, setup, and execution is in the ``sunbeamlib`` module. This module is located in ``src/sunbeamlib/`` with the root ``pyproject.toml`` configuring it. It has a number of different scripts, each located in its own file prefixed by ``script_``, and also has utility functions and classes for both the scripts and some portions of the pipeline that are particularly common.

workflow
========

The core of the work done by Sunbeam is handled by the Snakemake workflow, located in ``workflow/``. Once a project is setup properly by sunbeamlib, the workflow can be run with all the reproducibility benefits of Snakemake. The core of the workflow is defined in ``workflow/Snakefile``. Reference the Snakemake docs for help understanding Snakemake things better; they're very good. From this core Snakefile, we import more ``.smk`` files from ``workflow/rules/`` and ``extensions/sbx_*/``.

Important Variables
-------------------

Variables defined in the main Snakefile can be accessed throughout the workflow. Some important variables include:

- ``Samples``: Dict[str, Dict[str, str]] - A dictionary where keys are sample names and values are dictionaries of read pairs mapping to file paths (``Samples[sample] = {"1": r1, "2": r2}``).
- ``Pairs``: List[str] - Either ``["1", "2"]`` or ``["1"]`` depending on if the project is paired end or not.
- ``Cfg``: Dict[str, Dict[str, str]] - The YAML config converted into dictionary form.
- ``MIN_MEM_MB``: int - A minimum value of the number of megabytes of memory to request for each job. This will only apply for jobs that rely on Sunbeam to guess their memory requirements.
- ``MIN_RUNTIME``: int - A minimum value of the number of minutes to request for each job. This will only apply for jobs that rely on Sunbeam to guess their runtime requirements.
- ``HostGenomes``: List[str] - A list of host genomes that are used for decontaminating reads.
- ``HostGenomeFiles``: List[str] - A list of files with host genomes that are used for decontaminating reads (not to be confused with ``sbx_mapping``'s ``GenomeFiles`` variable, which it uses to track reference genome files).
- ``QC_FP``: Path - The Path to the project's quality control output directory.
- ``ASSEMBLY_FP``: Path - The Path to the project's assembly output directory.
- ``CLASSIFY_FP``: Path - The Path to the project's classification output directory.
- ``MAPPING_FP``: Path - The Path to the project's mapping output directory.
- ``BENCHMARK_FP``: Path - The Path to the project's benchmarking output directory.
- ``LOG_FP``: Path - The Path to the project's log output directory.

tests
=====

All tests are located in the ``tests/`` directory. The tests are run with pytest, and the tests are organized into subdirectories based on the module they are testing.

.github
=======

The ``.github/`` directory contains the configuration for GitHub Actions, which are used to run the tests on every push to the repository and manage releases. The configuration is in ``.github/workflows/``.