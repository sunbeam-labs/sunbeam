.. _install:

==========
install.sh
==========

Overview
========

This script enables users of sunbeam to easily install the necessary software and environments to run sunbeam.

Options
=======

All available options for the command line, used with ``./install.sh [options]``.

-e/--environment [arg]
+++++++++++++++++++++++++++++++
Environment to install to. Default: "sunbeam" followed by the version tag (e.g. sunbeam3.1.0). This version tag can get more complicated for non-release branches.

-s/--sunbeam_dir [arg]
+++++++++++++++++++++++++++++++
Location of Sunbeam source code. Default: root sunbeam directory. This should rarely be changed. The script is intended to be run from the root dir.

-c/--conda [arg]
+++++++++++++++++++++++++
Location of Conda installation. Default: $CONDA_PREFIX. If this variable doesn't point to your desired conda installation you can specify that here.

.. tip::
    The '--conda' argument should rarely be used. Having multiple conda installations or one where environment variables like $CONDA_PREFIX aren't set correctly is a recipe for trouble.

-u/--update [arg]
++++++++++++++++++++++++++
Update sunbeam [lib]rary, conda [env], or [all]. The sunbeam library is all code under the 'sunbeamlib/' directory. If you have major or incompatible changes to make to the environment you should consider creating a new one under a different name so that you always have a working version installed.

-v/--verbose
+++++++++++++++
Show subcommand output.

-d/--debug
+++++++++++++
Run in debug mode.

-h/--help
++++++++++++
Display help message.

-w/--version
+++++++++++++++
Display version information.