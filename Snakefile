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
import pandas
from io import StringIO
from pprint import pprint
from pathlib import Path

from snakemake.utils import update_config, listfiles
from snakemake.exceptions import WorkflowError

from sunbeamlib import build_sample_list
from sunbeamlib.config import *
from sunbeamlib.reports import *

def build_sample_from_barcode(bc_file):
    """"
    Build a list of samples from a barcode file
    :param bc_file: a Path to barcode file
    :returns: A dictionary of samples, with sample names as keys
    """
    with open(str(bc_file)) as f:
        lines = f.read().splitlines()
    ids = []
    for line in lines:
         ids.append(line.split("\t")[0])
    # todo: not sure about adding the path of actual reads
    Samples = dict((id,"paired") for id in ids)
    return Samples

if not config:
        raise SystemExit(
                "No config file specified. Run `sunbeam_init` to generate a "
                "config file, and specify with --configfile")

# ---- Substitute $HOME_DIR variable
#varsub(config)

# ---- Setting up config files and samples
Cfg = check_config(config)
Blastdbs = process_databases(Cfg['blastdbs'])

# ---- If the data_fp is empty, then read sample names from barcode_fp
if str(Cfg['all']['root']) != str(Cfg['all']['data_fp']):
     Samples = build_sample_list(Cfg['all']['data_fp'], Cfg['all']['filename_fmt'], Cfg['all']['exclude'])
else:
     Samples = build_sample_from_barcode(Cfg['all']['barcode_fp']) #todo: check existing of barcode_fp ?

# ---- Set up output paths for the various steps
QC_FP = output_subdir(Cfg, 'qc')
ASSEMBLY_FP = output_subdir(Cfg, 'assembly')
ANNOTATION_FP = output_subdir(Cfg, 'annotation')
CLASSIFY_FP = output_subdir(Cfg, 'classify')
MAPPING_FP = output_subdir(Cfg, 'mapping')

# ---- Change your workdir so .snakemake won't take all the space in the head node
workdir: str(Cfg['all']['output_fp'])

# ---- Targets rules
include: "rules/targets/targets.rules"


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
#include: "rules/classify/classify.rules"
#include: "rules/classify/clark.rules"
include: "rules/classify/kraken.rules"


# ---- Mapping rules
include: "rules/mapping/bowtie.rules"
include: "rules/mapping/kegg.rules"
#include: "rules/mapping/snap.rules"

# ---- Report rules
include: "rules/reports/reports.rules"


# ---- Rule all: run all targets
rule all:
    input: TARGET_ALL

rule samples:
    run:
        print("Samples found:")
        pprint(sorted(list(Samples.keys())))

onsuccess:
    print("Workflow finished, no error")
    shell("mail -s 'workflow finished' " + "zhaocy.dut@gmail.com < {log}")
onerror:
    print("An error occurred")
    shell("mail -s 'an error occurred' " + "zhaocy.dut@gmail.com < {log}")
