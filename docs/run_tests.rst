.. _run_tests:

==============
run_tests.bash
==============

.. contents::
   :depth: 2

Overview
========

This script manages the testing of sunbeam. It can be used to run a set of 
core tests that will verify a successful installation, run specific tests, or 
run all tests. The default behavior is to run the core set. It also gives the 
options to use a permanent directory for outputs (instead of a temporary one) 
or to use a pre-existing environment (instead of creating a temporary one).

Options
=======

All available options for the command line, used with ``tests/run_tests.sh [options]``.

-d [arg]
++++++++

Use [DIR] rather than a temporary directory (remains after tests finish).

-e [arg]
++++++++

Use a pre-existing Conda [ENV] rather than creating one (remains after tests finish).

-t [arg]
++++++++

Run a specific [TEST] from tests/test_suite.bash only or run [all] (default: runs core tests only).

-v
+++

Show subcommand output.

-h
+++

Display help message.

