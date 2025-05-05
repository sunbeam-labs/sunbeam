.. _commands:

================
Sunbeam Commands
================

Sunbeam
=======

.. argparse::
    :ref: sunbeam.scripts.sunbeam.main_parser
    :prog: sunbeam

Run
---

Usage examples
**************

1. To run all targets (not including extensions):
    ``sunbeam run --profile /path/to/project/``
2. To specify multiple targets:
    ``sunbeam run --profile /path/to/project/ all_decontam all_assembly all_annotation``
4. To run assembly on samples that have already been decontaminated:
    ``sunbeam run --profile /path/to/project/ --skip decontam all_assembly``

.. argparse::
    :ref: sunbeam.scripts.run.main_parser
    :prog: sunbeam run

Init
----

.. argparse::
    :ref: sunbeam.scripts.init.main_parser
    :prog: sunbeam init

Extend
------

.. argparse::
    :ref: sunbeam.scripts.extend.main_parser
    :prog: sunbeam extend

Config
------

.. argparse::
    :ref: sunbeam.scripts.config.main_parser
    :prog: sunbeam config
