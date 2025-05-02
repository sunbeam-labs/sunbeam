.. _commands:

================
Sunbeam Commands
================

Sunbeam
=======

.. argparse::
    :ref: sunbeam.scripts.sunbeam.main_parser
    :prog: sunbeam


.. argparse::
    :ref: sunbeam.scripts.run.main_parser
    :prog: sunbeam run
    :markdown:
    Run
    ===

    ## Usage examples

    1. To run all targets (not including extensions):
       ``sunbeam run --profile /path/to/project/``
    2. To specify multiple targets:
       ``sunbeam run --profile /path/to/project/ all_decontam all_assembly all_annotation``
    3. The equivalent of 2, using the deprecated ``--target_list`` option:
       ``sunbeam run --profile /path/to/project/ --target_list all_decontam all_assembly all_annotation``
    4. To run assembly on samples that have already been decontaminated:
       ``sunbeam run --profile /path/to/project/ --skip decontam all_assembly``


.. argparse::
    :ref: sunbeam.scripts.init.main_parser
    :prog: sunbeam init
    :markdown:
    Init
    ====
