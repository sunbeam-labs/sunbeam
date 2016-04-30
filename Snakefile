#
# Xenopeltis: an iridescent HTS pipeline
#
# Author: Erik Clarke <ecl@mail.med.upenn.edu>
# Created: 2016-04-28
#

import re
from pathlib import Path

from snakemake.utils import update_config, listfiles

Cfg_default = {
    'data_fp': '~/ext/100_SCID/103_Virome/data_files',
    'output_fp': '~/my_data/output',
    'filename_fmt': '{sample}_{rp}.fastq.gz',
    'subcores': 4,
    'genomes_fp': '~/ext/genomes',
    'bt2_index_fp': '~/ext/genomes/bowtie2_indices'
}

# These keys in the config file must point to valid paths
required_paths = ['data_fp', 'genomes_fp']

if Path('config.yaml').exists():
    configfile: 'config.yaml'
    update_config(Cfg_default, config)
else:
    print("No config.yaml found; using defaults")

# Deal with paths in the config file
Cfg = dict()
for k, v in Cfg_default.items():
    if k.endswith('_fp'):
        v = Path(v).expanduser()
        if k in required_paths and not v.exists():
            raise ValueError("Specified %s does not exist: %s" % (k, v))
        Cfg[k] = v
    else:
        Cfg[k] = v

# Build sample list
_all_files = list(listfiles(str(Cfg['data_fp']/Cfg['filename_fmt'])))
Samples = {t[1]['sample']: {} for t in _all_files}
for f in _all_files:
    rp = f[1]['rp'] if 'rp' in f[1].keys() else False
    if rp:
        Samples[f[1]['sample']][rp] = f[0]
        Samples[f[1]['sample']]['paired'] = True
    else:
        Samples[f[1]['sample']]['file'] = f[0]
        Samples[f[1]['sample']]['paired'] = False

## ==== Function definitions

include: "functions.snakefile"

## ==== Rules
        
rule all:
    input: expand(str(Cfg['output_fp']/'paired'/'{sample}.assembled.fastq'), sample = Samples.keys())

# ---- Quality control rules

include: "qc/qc.rules"

# ---- Host filtering rules

include: "qc/decontaminate.rules"

# ---- Assembly rules

include: "assembly/pairing.rules"
include: "assembly/assembly.rules"
    
# ---- Bowtie mapping rules

include: "bowtie.rules"


    
