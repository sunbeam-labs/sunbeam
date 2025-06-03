.. _faqs:

====
FAQs
====

A collection of common questions, issues, or points of confusion.

**I'm getting ``snakemake: error: argument --executor/-e: invalid choice: '_____' (choose from 'local', 'dryrun', 'touch')``. Why can't I use the ``--executor`` option?**

You're using the exectuor option properly, it's just that you haven't installed the executor plugin. Use ``pip`` to install it and you should be good to go (e.g. for Slurm ``pip install snakemake-executor-plugin-slurm``).

**I'm trying to use singularity but it keeps failing and complaining about running out of space. I know I have plenty of open disk space. Why is it running out?**

This is a known issue with singularity. It's not actually running out of space, it's just that the default location for the temporary directory is on a partition that is too small. You can change the location of the temporary directory by setting the ``SINGULARITY_TMPDIR`` and ``TMPDIR`` environment variables to a location with more space.

**A rule keeps failing with an error like "perl: error while loading shared libraries: libcrypt.so.1: cannot open shared object file: No such file or directory". What's going on?**

This is unfortunately a common issue with conda where shared libraries are either not installed or not properly loaded for packages that depend on them. There can be many causes and many fixes. You can start by searching the exact error message and seeing if there are any suggestions for how to solve it. Often it will involve installing the missing library with conda or installing the missing library with the system package manager. For example, the solution to the example error for me running sunbeam on a standard Amazon machine image (AMI) was to install the library using ``sudo yum install libxcrypt-compat`` in the snakemake-managed conda environment.

.. tip::
    
    If you can't find the answer to your problem here, try searching the `issue tracker <https://github.com/sunbeam-labs/sunbeam/issues>`_ on GitHub or posting a new issue.