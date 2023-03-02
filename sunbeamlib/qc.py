"""
Supporting functions for QC rules.
"""

import gzip
import re
from sunbeamlib.parse import parse_fastq, write_fastq


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
        for record in parse_fastq(f_in):
            trimmed_id = re.sub(suffix_pattern + "$", "", record[0])
            write_fastq((trimmed_id, record[1], record[2], record[3]), f_out)
    finally:
        f_in.close()
        f_out.close()


def filter_ids(fp_in, fp_out, ids, log):
    """Remove ids from FASTQ file.

    fp_in: path to input FASTQ
    fp_out: path to output FASTQ
    ids: list of ids to be removed
    """
    with gzip.open(fp_in, "rt") as f_in, gzip.open(fp_out, "wt") as f_out:
        for record in parse_fastq(f_in):
            if remove_pair_id(record[0]) not in ids:
                write_fastq(record, f_out)
            else:
                log.write(f"{record[0]} filtered\n")

def remove_pair_id(id, log):
    """Remove the 1 or 2 from a paired read ID
    
    id: id string
    """
    if id[:-2] == "/1" or id[:-2] == "/2":
        return id[:-2]
    space_split = id.split(" ")
    if len(space_split) == 2 and (space_split[1][0] == "1" or space_split[1][0] == "2"):
        return " ".join([space_split[0], space_split[1][1:]])
    log.write(f"Didn't find read pair ID in {id}")
    return id
