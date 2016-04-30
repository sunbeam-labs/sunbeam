#
# Xenopeltis: an iridescent HTS pipeline
#
# Author: Erik Clarke <ecl@mail.med.upenn.edu>
# Created: 2016-04-28
#

import re
from pathlib import Path
from snakemake.utils import update_config, listfiles

include: 'functions.snakefile'
configfile: 'config.yaml'

# Build config file, expanding and checking paths as required
Cfg = validate_paths(config)

QC_FP = Cfg['output_fp']/Cfg['qc_suffix']
ASSEMBLY_FP = Cfg['output_fp']/Cfg['assembly_suffix']
ANNOTATION_FP = Cfg['output_fp']/Cfg['annotation_suffix']

# Create sample list from parameters in config file
Samples = build_sample_list(Cfg['data_fp'], Cfg['filename_fmt'], Cfg['exclude'])

print(Samples.keys())

rule all:
    run:
        print("Xenopeltis: an iridescent HTS pipeline")
        print("For available commands, type `snakemake --list`")
        
# ---- Quality control rules
include: "qc/qc.rules"
include: "qc/decontaminate.rules"

# ---- Assembly rules
include: "assembly/pairing.rules"
include: "assembly/assembly.rules"

# ---- Viral annotation rules
include: "virome/annotation.rules"
    
# ---- Bowtie mapping rules
include: "bowtie.rules"
    

