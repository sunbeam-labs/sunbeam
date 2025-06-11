.. _extensions:

==================
Sunbeam Extensions
==================

Sunbeam is designed to be extensible. You can install and run extensions from our curated repository, other users' custom contributions, or your own. Extensions are simply Snakemake workflows that can be included on top of the core Sunbeam workflow. They can be used to add new functionality, integrate with other tools, or customize the pipeline to fit your needs.

.. tip::

    Extensions live in the ``extensions/`` directory under the Sunbeam root installation. You can also set the environment variable ``$SUNBEAM_EXTENSIONS`` to point to a different directory if you want to keep your extensions separate from the core Sunbeam installation.

This is now going to pick up where :ref:`quickstart` left off and show how to install and run ``sbx_kraken`` for taxonomic classification of your QCed and decontaminated reads.

Installation
************

To install the extension, run the following command:

.. code-block:: shell

   sunbeam extend sbx_kraken

Project setup
*************

You should already have a project set up but we need to update it with the new extension's configuration options. To do this, run the following command:

.. code-block:: shell

   sunbeam config update my_project/sunbeam_config.yml

In this case, this will add the following lines to your ``sunbeam_config.yml`` file:

.. code-block:: yaml

    sbx_kraken:
        kraken_db_fp: ''

You will need to edit this line to point to the Kraken database you want to use (find a Kraken database from their website or build your own).

Running the extension
*********************

Now that you have the extension installed and configured, you can run it with the following command:

.. code-block:: shell
   
   sunbeam run --profile my_project/ all_kraken

This will run the pipeline, recognizing that QC and decontam have already completed, and will run the Kraken extension on the decontaminated reads. The output will be stored in ``my_project/sunbeam_output/classify/kraken/``.