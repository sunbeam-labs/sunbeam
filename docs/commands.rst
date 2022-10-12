.. _commands:

================
Sunbeam Commands
================

.. contents::
   :depth: 2

sunbeam [-h/--help,-v/--version] <subcommand>

-h/--help: Display help.

-v/--version: Display version.

.. tip::

    This is the version given by `semantic_version` in sunbeamlib. For commits 
    and branches this may differ from the version given by `git` that is used 
    for things like environment naming.

init
====

sunbeam init [-h] [-f] [--output FILE] [--defaults FILE] [--template FILE] 
[--data_fp PATH] [--format STR] [--single_end] project_fp

-h: Display help.

-f: 

--output: 

--defaults: 

--template: 

--data_fp: 

--format: 

--single_end: 

project_fp: 

run
===

sunbeam run [-h] [-s PATH] -- <snakemake options>

-h: Display help.

-s: 

<snakemake options>: You can pass further arguments to Snakemake after `--`, 
e.g: `$ sunbeam run -- --cores 12`. See http://snakemake.readthedocs.io for 
more information.

config
======

sunbeam config [-h] {update,modify} ...

-h: Display help.

update
******

sunbeam config update [-h] [-t FILE] [--strict] [-i | -o FILE] config_file

-h: Display help.

-t: 

--strict: 

-i: 

-o: 

config_file: 

modify
******

sunbeam config modify [-h] [-s STR | -f FILE] [-i | -o FILE] config_file

-h: Display help.

-s: 

-f: 

-i: 

-o: 

config_file: 

list_samples 
============

sunbeam list_samples [-h] [-s] [-f STR] data_fp

-h: Display help.

-s: 

-f: 

data_fp: 

extend
======

sunbeam extend [-h] [-s PATH] github_url

-h: Display help.

-s/--sunbeam_dir: Path to sunbeam installation.

