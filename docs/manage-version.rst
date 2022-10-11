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

-l --list [arg]
---------

List all [installed] or all [available] versions of sunbeam.

.. _man-vm:
Manual Version Management
=========================



.. tip::

    This is the recommended method for developers of sunbeam. While the script 
    provides useful utilities, it can become a nuisance if you're making a lot 
    of updates to the code as it will keep trying to make new environments 
    when you switch.