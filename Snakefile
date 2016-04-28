#
# Xenopeltis: an iridescent HTS pipeline
#
# Author: Erik Clarke <ecl@mail.med.upenn.edu>
# Created: 2016-04-28
#

import re
from pathlib import Path
from pprint import pprint

from snakemake.utils import update_config, listfiles

Cfg_default = {
    'data_fp': '~/ext/100_SCID/103_Virome/data_files',
    'output_fp': '~/my_data/output',
    'filename_fmt': '{sample}_{rp}.fastq.gz',
    'subthreads': 4,
    'genomes_fp': '~/ext/genomes'
}

# These keys in the config file must point to valid paths
required_paths = ['data_fp', 'genomes_fp']

if Path('config.yaml').exists():
    configfile: 'config.yaml'
    Cfg_default.update_config(config)
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

## ==== Rule List ==== ##
        
rule all:
    input: expand(str(Cfg['output_fp']/'paired'/'{sample}.assembled.fastq'), sample = Samples.keys())

rule pair_reads:
    """Pairs reads with pear."""
    input:
        r1 = lambda wc: Samples[wc.sample]['R1'],
        r2 = lambda wc: Samples[wc.sample]['R2']
    output:
        str(Cfg['output_fp']/'paired'/'{sample}.assembled.fastq')
    params:
        out = str(Cfg['data_fp']/'paired'/'{sample}')
    shell:
        "pear -j {Cfg[subthreads]} -f {input.r1} -r {input.r2} -o {params.out}"
