import warnings
warnings.filterwarnings('ignore', '.*experimental.*')
from pathlib import Path
from collections import Counter, OrderedDict
import re
import os
import sys

import pandas
from io import StringIO
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord

from xml.etree.ElementTree import ParseError

def blast_summary(blast_files, blast_format="blast-xml"):
    """Summarize BLAST results from an set of BLAST output files."""
    for infile in blast_files:
        sample = Path(infile).stem
        try:
            for result in SearchIO.parse(infile, blast_format):
                if len(result.hits) > 0:
                    yield {
                        'sample': sample,
                        'query': result.id,
                        'hit':result.hits[0].id
                    }
        except ParseError:
            print("Skipping empty/malformed %s" % infile)
            continue

def blast_contig_summary(xml_files):
    return {r['query']: r['hit'] for r in blast_summary(xml_files)}

def blast_hits(blast_xml_files):
    return Counter(c['query'] for c in blast_summary(blast_xml_files))

def parse_trim_summary_paired(f):
    for line in f.readlines():
        if line.startswith('Input Read'):
            vals = re.findall('\D+\: (\d+)', line)
            keys = ('input', 'both_kept','fwd_only','rev_only','dropped')
            return(OrderedDict(zip(keys, vals)))

def parse_trim_summary_single(f):
    for line in f:
        if line.startswith('Input Read'):
            vals = re.findall('\D+\: (\d+)', line)
            keys = ('input', 'kept', 'dropped')
            return(OrderedDict(zip(keys, vals)))

def parse_decontam_log(f):
    keys = f.readline().rstrip().split('\t')
    vals = f.readline().rstrip().split('\t')
    return(OrderedDict(zip(keys,vals)))

def summarize_qual_decontam(tfile, dfile, paired_end):
    """Return a dataframe for summary information for trimmomatic and decontam rule"""
    tname = os.path.basename(tfile).split('.out')[0]
    dname = os.path.basename(dfile).split('.txt')[0]
    with open(tfile) as tf:
        with open(dfile) as jf:
            if paired_end:
                trim_data = parse_trim_summary_paired(tf)
            else:
                trim_data = parse_trim_summary_single(tf)
                
            decontam_data = parse_decontam_log(jf)
    sys.stderr.write("trim data: {}\n".format(trim_data))
    sys.stderr.write("decontam data: {}\n".format(decontam_data))
    return(pandas.DataFrame(OrderedDict(trim_data, **(decontam_data)), index=[tname]))

def parse_fastqc_quality(filename):
    with open(filename) as f:
        report = f.read()
    tableString = re.search(
        '\>\>Per base sequence quality.*?\n(.*?)\n\>\>END_MODULE',
        report, re.DOTALL).group(1)

    f_s = StringIO(tableString)
    df = pandas.read_csv(
        f_s, sep='\t', usecols=['#Base', 'Mean'], index_col='#Base')
    sample_name = os.path.basename(filename.split('_fastqc')[0])
    df.columns=[sample_name]
    f_s.close()
    return(df)

    


        
    
