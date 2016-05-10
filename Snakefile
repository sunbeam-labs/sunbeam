#
# Sunbeam: an iridescent HTS pipeline
#
# Author: Erik Clarke <ecl@mail.med.upenn.edu>
# Created: 2016-04-28
#

import re
import yaml
from pprint import pprint
from pathlib import Path

from snakemake.utils import update_config, listfiles

from sunbeam import build_sample_list
from sunbeam.config import *
from sunbeam.reports import *

configfile: 'config.yml'

# ---- Setting up config files and samples
Cfg = validate_paths(config)
blastdbs = process_databases(yaml.load(open('databases.yml')))
Samples = build_sample_list(Cfg['data_fp'], Cfg['filename_fmt'], Cfg['exclude'])

# ---- Rule all: show intro message
rule all:
    run:
        print("Samples found:")
        pprint(list(Samples.keys()))
        print("For available commands, type `snakemake --list`")

# ---- Quality control rules
include: "qc/qc.rules"
include: "qc/decontaminate.rules"


# ---- Assembly rules
include: "assembly/pairing.rules"
include: "assembly/assembly.rules"


# ---- Contig annotation rules
include: "annotation/annotation.rules"
include: "annotation/blast.rules"
include: "annotation/orf.rules"


# ---- Classifier rules
include: "classify/classify.rules"


# ---- Bowtie mapping rules
# include: "bowtie.rules"
