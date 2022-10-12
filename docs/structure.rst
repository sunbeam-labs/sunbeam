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
which then invokes the necessary snakemake commands. 

Sections
========

sunbeam/ (root directory)
-------------------------



docs/
-----



envs/
-----



extensions/
-----------



rules/
------



scripts/
--------



sunbeamlib/
-----------



tests/
------



