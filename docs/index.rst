.. Sunbeam documentation master file, created by
   sphinx-quickstart on Thu Jan 18 16:59:48 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.



==================
Welcome to Sunbeam
==================

.. image:: images/sunbeam_logo.gif
   :width: 120px
   :height: 120px
   :align: left

Sunbeam is a pipeline written in `snakemake <http://snakemake.readthedocs.io>`_
that simplifies and automates many of the steps in metagenomic sequencing
analysis. Sunbeam requires a reasonably modern GNU/Linux computer with bash, 
Python 3.1+, internet access (to retrieve dependencies), 4Gb of RAM, and at 
least 3Gb of disk space. RAM and disk space requirements may increase depending 
on the databases and tasks you choose to run, and the size of your data. For 
more information, check out the `Sunbeam paper in Microbiome 
<https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-019-0658-x>`_.

Sunbeam currently automates the following tasks:

* Quality control, including adaptor trimming, host read removal, and quality
  filtering;
* Decontamination of host-contaminated reads

Sunbeam was designed to be modular and extensible. We have 
:ref:`extensions` available that handle assembly, annotation, read alignments, taxonomic classifications, and more.

To get started, see our :ref:`quickstart`!

If you use Sunbeam in your research, please cite:

EL Clarke, LJ Taylor, C Zhao *et al.* Sunbeam: an
extensible pipeline for analyzing metagenomic
sequencing experiments. *Microbiome* 7:46 (2019)

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Contents:

   quickstart.rst
   usage.rst
   commands.rst
   structure.rst
   extensions.rst
   install.rst
   manage-version.rst
   run_tests.rst
   citation.rst

