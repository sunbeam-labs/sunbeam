"""
Supporting functions for QC rules.
"""

import gzip
from sunbeamlib.parse import parse_fastq, write_fastq


def filter_ids(fp_in, fp_out, ids, log):
    """Remove ids from FASTQ file.

    fp_in: path to input FASTQ
    fp_out: path to output FASTQ
    ids: list of ids to be removed
    """
    with gzip.open(fp_in, "rt") as f_in, open(fp_out, "w") as f_out:
        for record in parse_fastq(f_in):
            if any(id in record[0] for id in ids):
                log.write(f"{record[0]} filtered\n")
            else:
                write_fastq(record, f_out)


def remove_pair_id(id, log):
    """Remove the 1 or 2 from a paired read ID

    id: id string
    """
    id = id.strip()
    if id[-2:] == "/1" or id[-2:] == "/2":
        return id[:-2]
    
    # Assuming it's the newer id variant where komplexity removes the second half (containing pair number)
    return id
