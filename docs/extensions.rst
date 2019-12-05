.. _extensions:

==============
Sunbeam Extensions
==============

Sunbeam extensions allow you to add new features to the pipeline. You
might use an extension to carry out a new type of analysis, or to
produce a report from Sunbeam output files. Extensions can be
re-distributed and installed by other researchers to facilitate
reproducible analysis.

A Sunbeam extension consists of a directory that contains a
*rules* file. The *rules* file must have a name that
ends with the file extension ".rules".  The name of the extension is
derived from the name of this file.  It is customary for the name of
the extension to start with "sbx_", which is shorthand for "Sunbeam
extension."  Technically, only the *rules* file is needed for a
Sunbeam extension.  In practice, almost all extensions include other
files to facilitate the installation of software dependencies, to
specify parameters in the configuration file, and to give instructions
to users.

In Sunbeam version 3.0 and higher, extensions can be installed using the
``sunbeam extend`` command, followed by the GitHub URL of the 
extension you're installing. For example, to install an extension to
run the kaiju classifier, you would run::

    sunbeam extend https://github.com/sunbeam-labs/sbx_kaiju/

Sunbeam extensions are installed by placing the extension directory in
the ``extensions/`` subdirectory of the Sunbeam software.  Once the
extension is in place, Sunbeam will find the *rules* file and
incorporate the new functionality into the Sunbeam workflow.

Writing the rules file
======================

The rules file contains the code for one or more rules in Snakemake
format. These rules describe how to take files produced by Sunbeam and
do something more with them.  When Sunbeam is run, the Snakemake
workflow management system determines how to put the rules together in
the right order to produce the desired result.

Example: MetaSPAdes assembler
-----------------------------

We'll use the `sbx_metaspades_example
<https://github.com/sunbeam-labs/sbx_metaspades_example>`_ extension
to illustrate the basics.  This extension runs the MetaSPAdes assembly
program on each pair of host-filtered FASTQ files, producing a folder
with contigs for each sample.

The rules file contains two rules: one describes how to run the
MetaSPAdes program, and the other gives the full set of output folders
we expect to make.  Here is the rule that runs the program::

    rule run_metaspades:
        input:
            r1 = str(QC_FP/'decontam'/'{sample}_1.fastq.gz'),
            r2 = str(QC_FP/'decontam'/'{sample}_2.fastq.gz')
        output:
            str(ASSEMBLY_FP/'metaspades'/'{sample}')
        shell:
            "metaspades.py -1 {input.r1} -2 {input.r2} -o {output}"

The first line indicates a new rule named ``run_metaspades``.  The
``input`` and ``output`` sections contain patterns for file paths that
will be used by the command and produced after the command is run.
The command itself is given in the ``shell`` section.  Files from the
input and output sections are indicated in curly braces inside the
command.  In the case of multiple inputs, you can assign them
individual names, like ``r1`` and ``r2``, then refer to these names
inside the command.

In this example, the input and output filepaths are given as Python
code.  This is typical for rules in Sunbeam, because we are using info
from the user's configuration to determine the full filepath.  Let's
take ``str(QC_FP/'decontam'/'{sample}_1.fastq.gz')`` and see how it's
put together.  ``QC_FP`` is a Python variable, which gives the
absolute path to the directory containing quality control results in
Sunbeam.  This value is a special ``Path`` object in Python, which
means that we can add subdirectories using the division symbol ``/``.
As punishment for this convenience, we have to convert this ``Path``
object back to an ordinary string before it can be used in the rule.
The ``str()`` function accomplishes this.  Here, ``decontam`` is the
name of a subdirectory inside the main quality control directory. The
last part, ``{sample}_1.fastq.gz`` looks like the name of a gzipped
FASTQ file, but has something going on inside those curly braces.

The ``{sample}`` bit inside the input and output sections is the most
important part of this rule.  It's called a *wildcard*, and indicates
that for any sample, we expect to see certain files before and after
the command is run.  If a sample is named ``Abc``, Snakemake will look
for the files ``Abc_1.fastq.gz`` and ``Abc_2.fastq.gz`` before running
the command.  If these files don't exist, that's an error.  Likewise,
Snakemake will check for the subdirectory ``metaspades/Abc`` after the
command is run; if it's not there, Snakemake will stop and report the
discrepancy.  Wildcards allow you to write one rule that runs on all
the samples.

We have a standard format for files that can be used as input to your
extension.  Here are all the patterns you can use:

+-----------------------+----------------------------------------------------------------+
| Sequence data files   | Target                                                         |
+=======================+================================================================+
| Quality-controlled,   | str(QC_FP/'cleaned'/'{sample}_{rp}.fastq.gz')                  |
| non-decontaminated    | str(QC_FP/'cleaned'/'{sample}_1.fastq.gz')                     |
| sequences             | str(QC_FP/'cleaned'/'{sample}_2.fastq.gz')                     |
+-----------------------+----------------------------------------------------------------+
| Quality-controlled,   | str(QC_FP/'decontam'/'{sample}_{rp}.fastq.gz')                 |
| decontaminated        | str(QC_FP/'decontam'/'{sample}_1.fastq.gz')                    |
| sequences             | str(QC_FP/'decontam'/'{sample}_2.fastq.gz')                    |
+-----------------------+----------------------------------------------------------------+
| Contig sequences      | str(ASSEMBLY_FP/'contigs'/'{sample}-contigs.fa')               |
+-----------------------+----------------------------------------------------------------+
| Open reading frame    | str(ANNOTATION_FP/'genes'/'prodigal'/'{sample}_genes_nucl.fa') |
| nucleotide sequences  |                                                                |
+-----------------------+----------------------------------------------------------------+
| Open reading frame    | str(ANNOTATION_FP/'genes'/'prodigal'/'{sample}_genes_prot.fa') |
| protein sequences     |                                                                |
+-----------------------+----------------------------------------------------------------+

+-----------------------+-----------------------------------------------+
| Summary tables        | Target                                        |
+=======================+===============================================+
| Attrition from        | str(QC_FP/'reports'/'preprocess_summary.tsv') |
| decontamination and   |                                               |
| quality control       |                                               |
+-----------------------+-----------------------------------------------+
| Sequence              | str(QC_FP/'reports'/'fastqc_quality.tsv')     |
| quality scores        |                                               |
+-----------------------+-----------------------------------------------+
| Taxonomic assignments | str(CLASSIFY_FP/'kraken'/'all_samples.tsv')   |
| from Kraken           |                                               |
+-----------------------+-----------------------------------------------+

Very often, you will use one of these patterns directly for input to a
rule.  For output from a rule, you'll often use a pattern from the
tables above and modify it to suit your command.

Our example extension has one more rule, representing the final set of
files that should be produced after the extension is run.  This rule
has no output and no command; only input.  After Python evaluates the
code, the input to this rule is the full list of directories created
by MetaSPAdes for every sample::

    rule all_metaspades:
        input:
            expand(str(ASSEMBLY_FP/'metaspades'/'{sample}'),
                   sample=Samples.keys())

This rule is critical for the ``{sample}`` pattern to work inside the
Snakemake workflow management system.  To determine the names of the
samples, Snakemake *works backwards*, starting with the files you
*would like to produce* at the end of the workflow.  Snakemake does
not work forward; you can't give it a list of samples or assume that
it will match against input files already present.  This may seem
strange, but this way of working allows Snakemake to assemble a
workflow containing only the steps that are needed to make a
particular set of output files.

Fortunately, there is a basic pattern employed to write rules like
this.  Here, we take the output pattern from our other rule; this
gives the pattern for the files we'd like to have at the end.  Then,
we use a function called ``expand`` to generate the full list of
files.  The ``expand`` function expects to get a list of all possible
values for every wildcard in the filename.  Sunbeam provides two
variables for this purpose: ``Samples.keys()`` gives the full list of
sample names, and ``Pairs`` gives the values used for the forward and
reverse reads in the file.  Here, we give ``sample=Samples.keys()`` as
an additional argument to ``expand()``, and the function produces a
list of all the outputs we expect.

When the user runs the extension, they specify the rule name,
``all_metaspades``.  Using the full list of output directories,
Snakemake figures out what sample files it needs to use, figures out
what commands to run, runs the commands in parallel if possible, and
lets you know if there were any problems.

Example: a reproducible report
----------------------------

As another example, we'll look at an extension that takes standard
output from Sunbeam and produces a report.  The extension
`sbx_shallowshotgun_pilot
<https://github.com/junglee0713/sbx_shallowshotgun_pilot>`_ enables
researchers to re-run the analysis for a small methods comaprison.

To make a report from Sunbeam output files, the extension needs only
one rule.    ::

  rule make_shallowshotgun_report:
      input:
          kraken = str(CLASSIFY_FP/'kraken'/'all_samples.tsv'),
          preprocess = str(QC_FP/'preprocess_summary.tsv'),
          quality = str(QC_FP/'fastqc_quality.tsv'),
          sampleinfo = sunbeam_dir + '/extensions/sbx_shallowshotgun_pilot/data/sampleinfo.tsv'
      output:
          str(Cfg['all']['output_fp']/'reports/ShallowShotgun_Pilot_Report.html')
      script:
          'shallowshotgun_pilot_report.Rmd'

Here, the output is a single file path, and the path does not contain
any wildcards like ``{sample}``.  Therefore, Snakemake can work
backwards from the output file and figure out everything it needs; we
can use this rule as our final target when running Sunbeam.

The basic structure of the rule and most of the inputs should be
familiar from the previous example.  One of the inputs,
``sampleinfo``, does not come from Sunbeam, but is distributed with
the extension.  We know the filepath inside the extension is
``data/sampleinfo.tsv``, but we need to specify the entire path for
Snakemake to find the file.  To do this, we use the variable
``sunbeam_dir``, which points to the Sunbeam installation directory.
The extension must be located inside the ``extensions/`` subdirectory
to run.  From here, we know how to get to our file.  Because the value
of ``sunbeam_dir`` is an ordinary string, we use the ``+`` symbol to
add on the ``extensions/`` subdirectory, the directory name for the
extension, and the path to the file inside the extension directory.
This example shows how to refer to files inside the Sunbeam
installation directory.

In the output section, we need to specify a file path for the final
report.  Here, we use the configuration parameter
``Cfg['all']['output_fp']`` to get the base directory for output from
Sunbeam.  The value of this configuration parameter is a ``Path``
object, so we use the ``/`` symbol to add the rest of the filepath,
and surround the whole thing with the ``str()`` function.  Just as a
note, Snakemake will create the ``reports/`` subdirectory if needed,
so you don't have to worry about directories being present ahead of
time to accommodate your output files.

At the bottom of the rule, we write ``script`` instead of ``shell``,
because we'd like Snakemake to run a script instead of a shell
command.  Here, we give the name of a script in `R Markdown
<https://rmarkdown.rstudio.com/>`_ format.  The file path of the
script is given *relative to the rules file*, which is a little bit
different from all the other file paths in the rules file, but
convenient.

Inside the script, we need to access the input files given in the
rule.  Here is the part of the script that accesses the input file
paths and saves them as ordinary variables in R::

  sample_fp <- file.path(snakemake@input[["sampleinfo"]])
  preprocess_fp <- file.path(snakemake@input[["preprocess"]])
  quality_fp <- file.path(snakemake@input[["quality"]])
  kraken_fp <- file.path(snakemake@input[["kraken"]])

The `R Markdown tutorial
<https://rmarkdown.rstudio.com/lesson-1.html>`_ and `book
<https://bookdown.org/yihui/rmarkdown/>`_ are the best sources of
information on the report format, whereas the `R for data science book
<https://r4ds.had.co.nz/>`_ provides a good introduction to the R
programming languageas you might use it in the report.

Variables provided by Sunbeam
-----------------------------

Here is a table of all the Python variables provided by Sunbeam for
use in your extensions:

+-------------------+-------------+----------------------------------------------+
| Variable name     | Type        | Description                                  |
+-------------------+-------------+----------------------------------------------+
| ``QC_FP``         | Path        | Output directory for quality control files.  |
+-------------------+-------------+----------------------------------------------+
| ``ASSEMBLY_FP``   | Path        | Output directory for assembly files.         |
+-------------------+-------------+----------------------------------------------+
| ``ANNOTATION_FP`` | Path        | Output directory for gene annotation files.  |
+-------------------+-------------+----------------------------------------------+
| ``CLASSIFY_FP``   | Path        | Output directory for taxonomic               |
|                   |             | classification files.                        |
+-------------------+-------------+----------------------------------------------+
| ``Samples``       | Dictionary  | Key is the sample name, value is a dictionary|
|                   |             | with keys "1" and "2", values are the        |
|                   |             | the gzipped FASTQ files at the start of the  |
|                   |             | workflow.  For unpaired reads the value for  |
|                   |             | "2" is the empty string.                     |
+-------------------+-------------+----------------------------------------------+
| ``Pairs``         | List        | For paired reads, ["1", "2"]. For unpaired   |
|                   |             | reads, ["1"].                                |
+-------------------+-------------+----------------------------------------------+
| ``Cfg``           | Dictionary  | Parameters found in the configuration file.  |
|                   |             | For any parameter ending in "_fp", the value |
|                   |             | is converted to a Path object.  The most     |
|                   |             | commonly used parameter is                   |
|                   |             | ``Cfg['all']['output_dir']``, which gives the|
|                   |             | base output directory.                       |
+-------------------+-------------+----------------------------------------------+
| ``sunbeam_dir``   | String      | File path where Sunbeam is installed.        |
+-------------------+-------------+----------------------------------------------+

Further reading
---------------

We're only scratching the surface of what you can do with rules in
Snakemake.  The `official Snakemake documentation
<https://snakemake.readthedocs.io/en/stable/index.html>`_ gives
excellent instructions with more examples.

Software dependencies
=====================

If your extension requires additional software to be installed, you
can provide the names of `Conda packages <https://conda.io/docs/>`_
inside a file named ``requirements.txt``.  This file contains the
package names, one per line.  To install Conda packages in this file,
users of your extension will run the ``conda install`` command with
this file as an additonal argument::

    conda install --file requirements.txt

Configuration
=============

Your extension can include its own section in the configuration file.
To take advantage of this, you would write an example configuration
file named ``config.yml``. This file should contain only one
additional configuration section, specifying parameters for your
extension.  For example, the `sbx_coassembly
<https://github.com/sunbeam-labs/sbx_coassembly>`_ extension includes
two parameters: the number of threads to use, and the path to a file
with groups of samples to co-assemble.::

  sbx_coassembly:
    threads: 4
      group_file: ''

As of version 3.0, config options from extensions are automatically included
in config files made using ``sunbeam init`` and ``sunbeam config update``. This
functionality depends on the extension's configuration file being named
``config.yml``.

In version <3.0, users can copy this example section to the end of their
configuration file, using ``cat``::

  cat config.yml >> /path/to/user/sunbeam_config.yml

In your *rules* file, you can access parameters in the configuration
like this: ``Cfg['sbx_coassembly']['group_file']``.

The README file
===============

We recommend that you include a README file in your extension.  The
contents of the file should be in Markdown format, and the file should
be named ``README.md``.  Here's what you should cover in the README file:

1. A short summary of what your extension does
2. Any relevant citations
3. Instructions to install
4. Instructions to configure
5. Instructions to run

A good example to follow is the `sbx_coassembly
<https://github.com/sunbeam-labs/sbx_coassembly>`_ extension.

Publishing at sunbeam-labs.org
==============================

You are welcome to add your Sunbeam extensions to the directory at
`sunbeam-labs.org <https://sunbeam-labs.org/>`_.  To submit your
extension to the directory, please go to the `development page for
sunbeam-labs.org
<https://github.com/sunbeam-labs/sunbeam-labs.github.io>`_ and `open
an issue
<https://github.com/sunbeam-labs/sunbeam-labs.github.io/issues>`_ with
the GitHub URL of your extension. If you know Javascript, you can edit
the list at the top of the file ``main.js`` and `send us a pull request
<https://github.com/sunbeam-labs/sunbeam-labs.github.io/pulls>`_.
