.. _manage-version:

=================
manage-version.sh
=================

.. contents::
   :depth: 3

Overview
========

This script exists to enable users of sunbeam to easily switch between different 
versions of sunbeam; automatically managing the environments and git 
branches/tags.

.. tip::

    This script was implemented in v3.1.0, so using it to switch to any prior 
    versions will then require referring to the man-vm_ section to switch 
    versions again.

Options
=======

All available options for the command line, used with `./manage-version.sh [options]`.

-l [arg] OR --list [arg]
++++++++++++++++++++++++

List all [installed] or all [available] versions of sunbeam. The [installed] 
argument will search your local conda environments for the prefix 'sunbeam'. 
The [available] argument will list available release tags and developement 
branches.

-a OR --active
++++++++++++++

List environment for the code currently installed (active branch tag). I.e. 
this will get the current git tag and display the appropriate environment.

-c OR --clean
+++++++++++++

Remove all auxiliary sunbeam conda environments. These will typically be stored 
in `$SUNBEAM_DIR/.snakemake/` and can take up space, especially if you're a 
developer and make changes to environment files often.

.. tip::

    The next time running sunbeam after `./manage-version.sh --clean` will 
    take longer because it has to remake some of the cleaned environments.

-s [arg] OR --switch [arg]
++++++++++++++++++++++++++

Switch to a new version of sunbeam (install if not installed). This version 
argument can be 'dev', 'stable', any other branch name, or any version tag. 
A list of available versions can be listed with 
`./manage-version.sh -l available`.

-r [arg] OR --remove [arg]
++++++++++++++++++++++++++

Uninstall the specified version of sunbeam. A list of installed versions can 
be shown with `./manage-version.sh -l installed`.

-v OR --verbose
+++++++++++++++

Show subcommand output.

-d OR --debug
+++++++++++++

Run in debug mode.

-h OR --help
++++++++++++

Display help message.

.. _man-vm:
Manual Version Management
=========================



.. tip::

    This is the recommended method for developers of sunbeam. While the script 
    provides useful utilities, it can become a nuisance if you're making a lot 
    of updates to the code as it will keep trying to make new environments 
    when you switch.