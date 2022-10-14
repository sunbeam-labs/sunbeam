.. _install:

==========
install.sh
==========

.. contents::
   :depth: 2

Overview
========

This script enables users of sunbeam to easily install the necessary software 
and environments to run sunbeam. For the typical user, this script will be 
called to perform the initial install of sunbeam and thereafter any upgrades 
will be handled by the manage-version script. If you are doing development work 
on sunbeam you will likely make use of this script more often, in particular 
the --update argument.

Options
=======

All available options for the command line, used with ``./install.sh [options]``.

-e/--environment [arg]
+++++++++++++++++++++++++++++++

Environment to install to. Default: "sunbeam" followed by the version tag 
(e.g. sunbeam3.1.0). This version tag can get more complicated for non-release 
branches and can be shown with ``./manage-version.sh -a``.

-s/--sunbeam_dir [arg]
+++++++++++++++++++++++++++++++

Location of Sunbeam source code. Default: root sunbeam directory. This should 
rarely be changed. These scripts are intended to be run from the root dir.

-c/--conda [arg]
+++++++++++++++++++++++++

Location of Conda installation. Default: $CONDA_PREFIX. If this variable 
doesn't point to your desired conda installation you can specify that here.

.. tip::

    The '--conda' argument should rarely be used. Having multiple conda 
    installations or one where environment variables like $CONDA_PREFIX aren't 
    set correctly is a recipe for trouble.

-u/--update [arg]
++++++++++++++++++++++++++

Update sunbeam [lib]rary, conda [env], or [all]. The sunbeam library is all 
code under the 'sunbeamlib/' directory. If you have major or incompatible 
changes to make to the environment you should be creating a new one using 
``./manage-version.sh -s ENV_NAME``. This way you will maintain the option to 
easily switch back and forth if you make breaking changes in the new one.

-m/--no_mamba
++++++++++++++++

Don't use mamba in base environment as dependency solver. It is the default 
option to use mamba because it is considerably faster than conda in solving new 
environments. However it can also sometimes be a pain to install, especially 
with crowded 'base' environments.

-v/--verbose
+++++++++++++++

Show subcommand output.

-d/--debug
+++++++++++++

Run in debug mode.

-h/--help
++++++++++++

Display help message.