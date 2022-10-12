.. _structure:

==================
Software Structure
==================

.. contents::
   :depth: 2

Overview
========

Sunbeam is a snakemake pipeline with a python library acting as a wrapper. 
Calling `sunbeam [cmd] [args] [options]` is a call to this wrapper library 
which then invokes the necessary snakemake commands. The main Snakefile can be 
found in the root directory and it makes use of rules from `rules/` and 
`extensions/`, scripts from `scripts/`, and environments from `envs/`. Tests 
are managed by a script `tests/run_tests.bash` which collects test 
functions from `tests/test_suite.bash`. Documentation lives in `docs/` and is 
served by ReadTheDocs.

Sections
========

sunbeam/ (root directory)
-------------------------

The root sunbeam directory holds a few important files including 
`environment.yml`, `setup.py`, and `Snakefile`. The environment file defines 
the dependencies required to run sunbeam and is used to create the main sunbeam 
environment. The setup file defines the structure and dependencies of the 
sunbeamlib_ and makes it installable via pip. The snakefile manages all of the 
snakemake components of the pipeline.

.. tip::

    `environment.yml` defines the main sunbeam environment that you activate in 
    order to run the pipeline. Internally, sunbeam then manages a number of 
    other environments (defined in envs_) on a per-rule basis.

docs/
-----


.. _envs:
envs/
-----



extensions/
-----------



rules/
------



scripts/
--------



.. _sunbeamlib:
sunbeamlib/
-----------



tests/
------



Hidden Directories
------------------

.circleci/
++++++++++



.github/
++++++++



.snakemake/
+++++++++++

This directory is created the first time you run sunbeam. 