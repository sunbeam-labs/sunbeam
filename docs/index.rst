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
Python 2.6+, internet access (to retrieve dependencies), 4Gb of RAM, and at 
least 3Gb of disk space. RAM and disk space requirements may increase depending 
on the databases and tasks you choose to run, and the size of your data.

Sunbeam currently automates the following tasks:

* Quality control, including adaptor trimming, host read removal, and quality
  filtering;
* Taxonomic assignment of reads to databases using `Kraken
  <https://github.com/DerrickWood/kraken>`_;
* Assembly of reads into contigs using `Megahit
  <https://github.com/voutcn/megahit>`_;
* Contig annotation using BLAST[n/p/x];
* Mapping of reads to target genomes; and
* ORF prediction using `Prodigal <https://github.com/hyattpd/Prodigal>`_

Sunbeam was designed to be modular and extensible. We have a few pre-built
:ref:`extensions` available that handle visualization tasks, including contig
assembly graphs, read alignments, and taxonomic classifications.

To get started, see our :ref:`quickstart`!

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Contents:

   quickstart.rst
   usage.rst
   extensions.rst

