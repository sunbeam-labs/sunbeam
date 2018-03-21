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
from pathlib import Path, PurePath

from snakemake.utils import update_config, listfiles
from snakemake.exceptions import WorkflowError

from sunbeamlib import build_sample_list, read_seq_ids
from sunbeamlib.config import *
from sunbeamlib.reports import *

# Load config file
if not config:
    raise SystemExit(
        "No config file specified. Run `sunbeam init` to generate a "
        "config file, and specify with --configfile")

# Check for major version compatibility
pkg_major, cfg_major = check_compatibility(config)
if pkg_major > cfg_major:
    raise SystemExit(
        "\nThis config file was created with an older version of Sunbeam"
        " and may not be compatible. Create a new config file using"
        "`sunbeam init` and update it using `sunbeam_mod_config`\n")
elif pkg_major < cfg_major:
    raise SystemExit(
        "\nThis config file was created with an newer version of Sunbeam"
        " and may not be compatible. Create a new config file using "
        "`sunbeam init` and update it using `sunbeam_mod_config`\n")

# Load extensions
sbxs = list(listfiles("extensions/{sbx_folder}/{sbx}.rules"))
for sbx in sbxs:
    sys.stderr.write("Found extension {sbx} in folder {sbx_folder}\n".format(**sbx[1]))

# Setting up config files and samples
Cfg = check_config(config)
Blastdbs = process_databases(Cfg['blastdbs'])
Samples = build_sample_list(Cfg['all']['samplelist_fp'], Cfg['all']['paired_end'])
print(Samples)

# Collect host (contaminant) genomes
sys.stderr.write("Collecting host/contaminant genomes... ")
HostGenomeFiles = [f for f in Cfg['qc']['host_fp'].glob('*.fasta')]
if not HostGenomeFiles:
    sys.stderr.write(
        "\n\nWARNING: No files detected in host genomes folder ({}). "
        "If this is not intentional, make sure all files end in "
        ".fasta and the folder is specified correctly.\n\n".format(
            Cfg['qc']['host_fp']
        ))
HostGenomes = {Path(g.name).stem: read_seq_ids(Cfg['qc']['host_fp'] / g) for g in HostGenomeFiles}
sys.stderr.write("done.\n")

sys.stderr.write("Collecting target genomes... ")
GenomeFiles = [f for f in Cfg['mapping']['genomes_fp'].glob('*.fasta')]
GenomeSegments = {PurePath(g.name).stem: read_seq_ids(Cfg['mapping']['genomes_fp'] / g) for g in GenomeFiles}
sys.stderr.write("done.\n")

# ---- Change your workdir to output_fp
workdir: str(Cfg['all']['output_fp'])

# ---- Set up output paths for the various steps
QC_FP = output_subdir(Cfg, 'qc')
ASSEMBLY_FP = output_subdir(Cfg, 'assembly')
ANNOTATION_FP = output_subdir(Cfg, 'annotation')
CLASSIFY_FP = output_subdir(Cfg, 'classify')
MAPPING_FP = output_subdir(Cfg, 'mapping')


# ---- Targets rules
include: "rules/targets/targets.rules"


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
include: "rules/classify/kraken.rules"


# ---- Mapping rules
include: "rules/mapping/mapping.rules"


# ---- Reports rules
include: "rules/reports/reports.rules"

for sbx_path, wildcards in sbxs:
    include: sbx_path
        

# ---- Rule all: run all targets
rule all:
    input: TARGET_ALL

rule samples:
    message: "Samples to be processed:"
    run:
        [print(sample) for sample in sorted(list(Samples.keys()))]


