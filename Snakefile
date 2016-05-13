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
Cfg = check_config(config)
print(Cfg['all']['data_fp'].exists())
Blastdbs = process_databases(yaml.load(open('databases.yml')))
Samples = build_sample_list(Cfg['all']['data_fp'], Cfg['all']['filename_fmt'], Cfg['all']['exclude'])

# ---- Set up output paths for the various steps
QC_FP = output_subdir(Cfg, 'qc')
ASSEMBLY_FP = output_subdir(Cfg, 'assembly')
ANNOTATION_FP = output_subdir(Cfg, 'annotation')
CLASSIFY_FP = output_subdir(Cfg, 'classify')

# ---- Rule all: show intro message
rule all:
    run:
        print("Samples found:")
        pprint(list(Samples.keys()))
        print("For available commands, type `snakemake --list`")

# ---- Quality control rules
include: "rules/qc/qc.rules"
include: "rules/qc/decontaminate.rules"


# ---- Assembly rules
include: "rules/assembly/assembly.rules"


# ---- Contig annotation rules
include: "rules/annotation/annotation.rules"
include: "rules/annotation/blast.rules"
include: "rules/annotation/orf.rules"


# ---- Classifier rules
include: "rules/classify/classify.rules"
