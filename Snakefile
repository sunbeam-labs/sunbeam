#
# Sunbeam: an iridescent HTS pipeline
#
# Author: Erik Clarke <ecl@mail.med.upenn.edu>
# Created: 2016-04-28
#

import re
import sys
import yaml
import configparser
from pprint import pprint
from pathlib import Path

from snakemake.utils import update_config, listfiles
from snakemake.exceptions import WorkflowError

from sunbeam import build_sample_list
from sunbeam.config import *
from sunbeam.reports import *

#if not config:
#    raise SystemExit("\nNo config file specified; specify a config file using --configfile\n")

configfile: "example_config.yml"

# ---- Setting up config files and samples
Cfg = check_config(config)
print(Cfg['assembly'])
Blastdbs = process_databases(Cfg['blastdbs'])
Samples = build_sample_list(Cfg['all']['data_fp'], Cfg['all']['filename_fmt'], Cfg['all']['exclude'])

# ---- Set up output paths for the various steps
QC_FP = output_subdir(Cfg, 'qc')
ASSEMBLY_FP = output_subdir(Cfg, 'assembly')
ANNOTATION_FP = output_subdir(Cfg, 'annotation')
CLASSIFY_FP = output_subdir(Cfg, 'classify')
MAPPING_FP = output_subdir(Cfg, 'mapping')

# ---- Targets rules
include: "rules/targets/targets.rules"

# ---- Rule all: show intro message
rule all:
    input: TARGET_FP
    run:
        print("Samples found:")
        pprint(sorted(list(Samples.keys())))

# ---- Quality control rules
include: "rules/qc/qc.rules"
include: "rules/qc/decontaminate.rules"


# ---- Assembly rules
include: "rules/assembly/assembly.rules"
include: "rules/assembly/pairing.rules"


# ---- Antibiotic resistance gene rules
include: "rules/abx/abx_genes.rules"

# ---- Contig annotation rules
include: "rules/annotation/annotation.rules"
include: "rules/annotation/blast.rules"
include: "rules/annotation/orf.rules"


# ---- Classifier rules
include: "rules/classify/classify.rules"
include: "rules/classify/clark.rules"
include: "rules/classify/kraken.rules"


# ---- Mapping rules
include: "rules/mapping/bowtie.rules"
#include: "rules/mapping/snap.rules"
