.. _commands:

================
Sunbeam Commands
================

.. contents::
   :depth: 2

.. code-block:: shell
    
    sunbeam [-h | -v] <subcommand>

.. code-block:: shell

    -h/--help: Display help.
    -v/--version: Display version.

.. tip::

    This is the version given by ``semantic_version`` in sunbeamlib. For commits 
    and branches this may differ from the version given by ``git`` that is used 
    for things like environment naming.

init
====

Initialize a new Sunbeam project in a given directory, creating a new config 
file and (optionally) a sample list.

.. code-block:: shell

    sunbeam init [-h] [-f] [--output FILE] [--defaults FILE] [--template FILE] [--data_fp PATH] [--format STR] [--single_end] project_fp

.. code-block:: shell

    -h/--help: Display help.
    -f/--force: Overwrite files if they already exist.
    --output: Name of config file (default: sunbeam_config.yml).

Config file options:

.. code-block:: shell

    --defaults: Set of default values to use to populate config file.
    --template: Custom config file template, in YAML format.

Sample list options:

.. code-block:: shell

    --data_fp: Path to folder containing ``.fastq.gz`` files.
    --format: Filename format for --data_fp (default: guessed).
    --single_end: Fastq files are in single-end, not paired-end, format for --data_fp.
    project_fp: Project directory (will be created if it does not exist).

run
===

Executes the Sunbeam pipeline by calling Snakemake.

.. code-block:: shell

    sunbeam run [-h] [-s PATH] -- <snakemake options>

.. code-block:: shell

    -h/--help: Display help.
    -s/--sunbeam_dir: Path to sunbeam installation.
    <snakemake options>: You can pass further arguments to Snakemake after ``--``, e.g: ``$ sunbeam run -- --cores 12``. See http://snakemake.readthedocs.io for more information.

config
======

.. code-block:: shell

    sunbeam config [-h] {update,modify} ...

.. code-block:: shell

    -h/--help: Display help.

update
******

Updates a config file to be compatible with the active version of sunbeam.

.. code-block:: shell

    sunbeam config update [-h] [-t FILE] [--strict] [-i | -o FILE] config_file

Usage examples:

1. To update a config file in place:
    ``sunbeam config update -i my_config.yml``
2. To write an update copy to a new file:
    ``sunbeam config update old_config.yml -o new_config.yml``

.. code-block:: shell

    -h/--help: Display help.
    -t/--template: Path to custom config file template, in YAML format.
    --strict: Remove keys that no longer exist in the new config file.
    -i/--in_place: Alters config file in place.
    -o/--out: Where to write modified config file.
    config_file: Existing config file to update.

modify
******

Modifies a config file with the specified changes.

.. code-block:: shell

    sunbeam config modify [-h] [-s STR | -f FILE] [-i | -o FILE] config_file

Usage examples:

1. To apply a set of defaults to an existing config file in place:
    ``sunbeam config modify -i -f defaults.yml my_config.yml``
2. To change a single key:value pair in the 'mapping' section:
    ``sunbeam config modify -i -s 'mapping: {keep_unaligned: True}'``

.. code-block:: shell

    -h: Display help.
    -s/--str: YAML string (e.g. 'blast: {threads: 4}').
    -f/--file: YAML file with new config values.
    -i/--in_place: Alters config file in place.
    -o/-out: Where to write modified config file.
    config_file: Existing config file to modify.

list_samples 
============

List the samples found in the specified directory.

.. code-block:: shell

    sunbeam list_samples [-h] [-s] [-f STR] data_fp

.. code-block:: shell

    -h/---help: Display help.
    -s/--single_end: Samples are single-end (not paired-end).
    -f/--format: Filename format (e.g. {sample}_R{rp}.fastq.gz) (default: guessed).
    data_fp: Path to folder containing reads.

extend
======

Install the extension at the given URL.

.. code-block:: shell

    sunbeam extend [-h] [-s PATH] github_url

.. code-block:: shell

    -h/--help: Display help.
    -s/--sunbeam_dir: Path to sunbeam installation.
