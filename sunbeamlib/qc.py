"""
Supporting functions for QC rules.
"""

import gzip
import re
from Bio import SeqIO

def strip_seq_id_suffix(fp_in, fp_out, suffix_pattern):
    """Remove sequence ID suffix from entries in a FASTQ file.

    fp_in: path to input FASTQ
    fp_out: path to output FASTQ
    suffix_pattern: regular expression for suffix to remove, e.g., "/[12]".
    This will be anchored to the end of the ID string.
    """
    try:
        if fp_in.endswith(".gz"):
            f_in = gzip.open(fp_in, "rt")
        else:
            f_in = open(fp_in, "rU")
        if fp_out.endswith(".gz"):
            f_out = gzip.open(fp_out, "wt")
        else:
            f_out = open(fp_out, "w")
        for record in SeqIO.parse(f_in, "fastq"):
            record.id = re.sub(suffix_pattern + "$", "", record.id)
            SeqIO.write(record, f_out, "fastq")
    finally:
        f_in.close()
        f_out.close()
