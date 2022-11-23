#
# Sunbeam: an iridescent HTS pipeline
#
# Author: Erik Clarke <ecl@mail.med.upenn.edu>
# Created: 2016-04-28
#
import os
import re
import sys
import yaml
import configparser

from pprint import pprint
from pathlib import Path, PurePath

from snakemake.utils import update_config, listfiles
from snakemake.exceptions import WorkflowError

from sunbeamlib import load_sample_list, read_seq_ids
from sunbeamlib.config import *
from sunbeamlib.reports import *

# Disallow slashes in our sample names during Snakemake's wildcard evaluation.
# Slashes should always be interpreted as directory separators.
wildcard_constraints:
  sample="[^/]+"

# Load config file
if not config:
    raise SystemExit(
        "No config file specified. Run `sunbeam init` to generate a "
        "config file, and specify with --configfile")

sunbeam_dir = ""
try:
    sunbeam_dir = os.environ["SUNBEAM_DIR"]
except KeyError:
    raise SystemExit(
        "$SUNBEAM_DIR environment variable not defined. Are you sure you're "
        "running this from the Sunbeam conda env?")

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
sbxs = list(listfiles(sunbeam_dir+"/extensions/{sbx_folder}/{sbx}.rules")) + list(listfiles(sunbeam_dir+"/extensions/{sbx_folder}/{sbx}.smk")) #commented to break test
for sbx in sbxs:
    sys.stderr.write("Found extension {sbx} in folder {sbx_folder}\n".format(**sbx[1]))

# Setting up config files and samples
Cfg = check_config(config)
Blastdbs = process_databases(Cfg['blastdbs'])
Samples = load_sample_list(Cfg['all']['samplelist_fp'], Cfg['all']['paired_end'], Cfg['all']['download_reads'], Cfg["all"]['root']/Cfg['all']['output_fp'])
Pairs = ['1', '2'] if Cfg['all']['paired_end'] else ['1']


# Collect host (contaminant) genomes
sys.stderr.write("Collecting host/contaminant genomes... ")
if Cfg['qc']['host_fp'] == Cfg['all']['root']:
    HostGenomeFiles = []
else:
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
if Cfg['mapping']['genomes_fp'] == Cfg['all']['root']:
    GenomeFiles = []
    GenomeSegments = {}
else:
    GenomeFiles = [f for f in Cfg['mapping']['genomes_fp'].glob('*.fasta')]
    GenomeSegments = {PurePath(g.name).stem: read_seq_ids(Cfg['mapping']['genomes_fp'] / g) for g in GenomeFiles}
sys.stderr.write("done.\n")

# ---- Change your workdir to output_fp
workdir: str(Cfg['all']['output_fp'])

# ---- Set up output paths for the various steps
DOWNLOAD_FP = output_subdir(Cfg, 'download')
QC_FP = output_subdir(Cfg, 'qc')
ASSEMBLY_FP = output_subdir(Cfg, 'assembly')
ANNOTATION_FP = output_subdir(Cfg, 'annotation')
CLASSIFY_FP = output_subdir(Cfg, 'classify')
MAPPING_FP = output_subdir(Cfg, 'mapping')


# ---- Targets rules
include: "rules/targets/targets.smk"


# ---- Quality control rules
include: "rules/qc/qc.smk"
include: "rules/qc/decontaminate.smk"


# ---- Assembly rules
include: "rules/assembly/assembly.smk"


# ---- Contig annotation rules
include: "rules/annotation/orf.smk"


for sbx_path, wildcards in sbxs:
    include: sbx_path
        

# ---- Rule all: run all targets
rule all:
    input: TARGET_ALL

rule samples:
    message: "Samples to be processed:"
    run:
        [print(sample) for sample in sorted(list(Samples.keys()))]
