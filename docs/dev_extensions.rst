.. _dev_extensions:

=============================
Developing Sunbeam Extensions
=============================

Sunbeam extensions allow you to add new features to the pipeline. You might use an extension to carry out a new type of analysis, or to produce a report from Sunbeam output files. Extensions can be re-distributed and installed by other researchers to facilitate reproducible analysis.

This page will go deep into the weeds on how to develop extensions for Sunbeam. If you're looking for how to install/run them, check out :ref:`extensions`. Make sure you've read the :ref:`dev` page first so that you have a good grounding in the internals of Sunbeam. And the more you're familiar with Snakemake and Python, the better (both have great documentation).

Starting a new extension
========================

A Sunbeam extension consists of a directory that contains a rules file. The rules file must have a name that ends with the file extension ".smk". It is customary for the name of the extension to start with "sbx_", which is shorthand for "Sunbeam extension."  Technically, only the rules file is needed for a Sunbeam extension. In practice, almost all extensions include other files to facilitate the installation of software dependencies, to specify parameters in the configuration file, and to give instructions to users.

Using ``sbx_template``
----------------------

We've created a GitHub template for creating new Sunbeam extensions easily. To use the template, go to the `sbx_template <https://github.com/sunbeam-labs/sbx_template>`_ repo and click the Use this template button in the top right. Once you've created the new repository, wait a minute for the CI to update all the names and links, then clone the repository to your computer. You can then edit the files in the repository to create your extension.

Writing the rules file
----------------------

The rules file contains the code for one or more rules in Snakemake format. These rules describe how to take files produced by Sunbeam and do something more with them.  When Sunbeam is run, the Snakemake workflow management system determines how to put the rules together in the right order to produce the desired result.

Example: An extension for AwesomeTool
-------------------------------------

Say we have a tool we want to include in Sunbeam called AwesomeTool, that provides strain level taxonomic classification of every read and then provides a human-readable report of how interesting each one is in descending order of interest. (unfortunately, AwesomeTool doesn't exist outside of this example). AwesomeTool has a command line interface (CLI) that takes gzipped fastq files and a database file as input, and produces a gzipped output file. The command looks like:

.. code-block:: shell
    AwesomeTool -i {input} -d {database} -o {output}

We'll start our development effort by updating the rules file (this assumes you're using ``sbx_template``; if not, follow along with those files in mind). We want to remove the dummy rules from the template ``sbx_awesome_tool.smk`` and replace them with our own:

.. code-block:: python
    AWESOME_FP = CLASSIFY_FP / 'awesome'

    rule all_awesome_tool:
        input:
            expand(AWESOME_FP / "{sample}.awesome", sample=Samples.keys())

    rule run_awesome_tool:
        input:
            r1 = QC_FP / "decontam" / "{sample}_1.fastq.gz",
            r2 = QC_FP / "decontam" / "{sample}_2.fastq.gz",
            db = Cfg['all']['awesome_db']
        output:
            AWESOME_FP / "{sample}.awesome"
        conda:
            "awesome_env.yaml"
        shell:
            """
            AwesomeTool -i {input.r1} -d {input.db} -o {output}
            """

This gives us our target (``all_awesome_tool``) and the rule that runs the program (``run_awesome_tool``). The target rule is a list of all the output files we expect to produce. The runner rule describes how to run the program, including the input files, output files, and command to run. You might notice that there are a few things here we will need to define to supplement these rules. For starters, we need to define the ``awesome_db`` option in the config. To do this, we will modify the ``config.yml`` file:

.. code-block:: yaml
    sbx_awesome_tool:
        awesome_db: /path/to/awesome_db

We also need to specify the environment that AwesomeTool is installed in. To do this, we will modify ``envs/sbx_awesome_tool_env.yml``:

.. code-block:: yaml
    name: awesome_env
    channels:
        - bioconda
        - conda-forge
    dependencies:
        - awesome_tool

And with that, we have a working extension! You can run it with the command:

.. code-block:: shell
    sunbeam run --profile /path/to/project/ all_awesome_tool

We can go check out our output files in ``/path/to/project/sunbeam_output/classify/awesome``.

Example: A reproducible report
------------------------------

After running ``sbx_awesome_tool`` a while, you might realize you're compiling the same standard report over and over again. To address this, we can create a new rule that produces a summary report on the results of AwesomeTool on each sample. To do this, we will add a new rule to our rules file and adjust our target rule:

.. code-block:: python
    rule all_awesome_tool:
        input:
            AWESOME_FP / "awesome_report.tsv",

    ...

    rule make_awesome_report:
        input:
            awesome = expand(AWESOME_FP / "{sample}.awesome", sample=Samples.keys())
        output:
            report = AWESOME_FP / "awesome_report.tsv"
        script:
            "scripts/awesome_report.R"

We choose R because we're familiar with it (but you should choose what you're familiar with). The ``script`` command points to a file ``scripts/awesome_report.R`` that we now have to make:

.. code-block:: r
    awesome_files <- snakemake@input$awesome
    awesome_report <- snakemake@output$report

    # Load the awesome files and make a report
    awesome_data <- lapply(awesome_files, read.table, header=TRUE)
    awesome_report_data <- do.call(rbind, awesome_data)
    awesome_report_data <- awesome_report_data[order(awesome_report_data$interest, decreasing=TRUE), ]
    awesome_report_data <- awesome_report_data[1:10, ]

    # Write the report
    write.table(awesome_report_data, awesome_report, sep="\t", row.names=FALSE, col.names=TRUE)

Note how the snakemake rule variables are accessed through the injected ``snakemake`` object. Access patterns differ accross languages but this is the general idea. After running the pipeline, we can check out our report in ``/path/to/project/sunbeam_output/classify/awesome/awesome_report.html``.

Example: Complicated use case
-----------------------------

AwesomeTool is great, but now we realize the standard CLI doesn't offer enough flexibility for what we want to do. We'll have to make a custom script to run it, importing some of the internals from the AwesomeTool package. To do this, we will need to modify our runner rule:

.. code-block:: python
    rule run_awesome_tool:
        input:
            r1 = QC_FP / "decontam" / "{sample}_1.fastq.gz",
            r2 = QC_FP / "decontam" / "{sample}_2.fastq.gz",
            db = Cfg['all']['awesome_db']
        output:
            AWESOME_FP / "{sample}.awesome"
        conda:
            "awesome_env.yaml"
        script:
            "scripts/run_awesome_tool.py"

And now we make the script:

.. code-block:: python
    import awesome_tool

    def run_awesome_tool(r1, r2, db):
        data = awesome_tool.load_data(r1, r2)
        results = awesome_tool.run_analysis(data, db)
        awesome_tool.save_results(results, snakemake.output[0])

    run_awesome_tool(snakemake.input.r1, snakemake.input.r2, snakemake.input.db)

Testing the extension
---------------------

Testing is a tricky thing with Snakemake. Because it is a workflow language, most of the code we write is just putting pieces together. So if we took the orthodox functional programming approach, we would test that the pipeline compiles (i.e. that we can run a dryrun without errors) and then claim that it will work because we can trust all the pieces to be tested themselves. So often we will *start* our testing with a dryrun:

.. code-block:: python
    import pytest
    import subprocess as sp
    from pathlib import Path

    @pytest.fixture
    def project(tmp_path):
        project_fp = tmp_path / "project"
        reads_fp = Path(".tests/data/reads").resolve()
        db_fp = Path(".tests/data/db").resolve()

        sp.check_output(["sunbeam", "init", "--data_fp", str(reads_fp), str(project_fp)])
        sp.check_output(["sunbeam", "config", "--modify", "sbx_awesome_tool: {awesome_db: " + str(db_fp) + "}"])

        return project_fp

    def test_dryrun(project):
        # Run a dryrun of the pipeline
        sp.check_output(["sunbeam", "run", "--profile", str(project), "--directory", str(project.parent), "--dryrun", "all_awesome_tool"])

.. tip::
    Adding the ``--directory`` flag pointing to a temp directory is a good idea for running these tests in a non-CI environment. This will prevent a bunch of snakemake garbage from being left behind in your working directory and make the tests more reliable.

This is nice, but in practice, we dont' really trust the components of our pipeline enough to stop here. If we have the ability to create a small enough test case to run the whole extension, that is an ideal case for testing everything easily. For example, here we can take our small set of synthetic reads in ``.tests/data/reads`` and a small dummy database we created with AwesomeTool's awesome database building utiltity in ``.tests/data/db``. We can then run the pipeline on these small files and check that the output is what we expect. This is a bit of a pain to set up, but it is worth it in the end. We can do this by adding this to our test file:

.. code-block:: python
    def test_run(project):
        sp.check_output(["sunbeam", "run", "--profile", str(project), "--directory", str(project.parent), "all_awesome_tool"])

        report_fp = project / "sunbeam_output" / "classify" / "awesome" / "awesome_report.tsv"

        assert report_fp.exists()
        assert report_fp.stat().st_size > 0
        with open(report_fp) as f:
            header = f.readline()
            first_line = f.readline()
            interest_index = header.split("\t").index("interest")
            assert float(first_line.split("\t")[interest_index]) > 0.5

If you don't have the ability to run all the rules of your extension in a CI environment (some tools just won't work on that small a dataset), you'll have to get more creative. You can test individual rules, explore Snakemake's unit testing framework, or test script logic directly.

.. tip::
    If you have enough custom logic that it needs testing, consider moving it to an extension library. This will allow you to test the logic separately from Snakemake's machinations. For instance, make a Python function for filtering reads if they have less than 20 A's. We would start by writing the function in the extension library, then calling it from the Snakemake script file (passing the snakemake input, output, etc objects), then testing the function.

The extension comes with a ``.github`` directory that defines CI workflows for automatically running tests.

Releasing the extension
-----------------------

Releases for Sunbeam extensions are pretty mellow. We don't have a proper repository or enforced release structure. You could never release your extension and it would still work fine. But if you want to maintain good versioning practices and developer hygiene, you can use GitHub's release system and update the version number in ``VERSION``.