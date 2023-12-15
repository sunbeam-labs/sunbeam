"""
Supporting functions for QC rules.
"""

import gzip
import re
from pathlib import Path
from sunbeamlib.parse import parse_fastq, write_fastq
from typing import List, TextIO


def filter_ids(fp_in: Path, fp_out: Path, ids: List[str], log: TextIO) -> None:
    """Remove ids from FASTQ file.

    fp_in: path to input FASTQ
    fp_out: path to output FASTQ
    ids: list of ids to be removed
    """
    with gzip.open(fp_in, "rt") as f_in, gzip.open(fp_out, "wt") as f_out:
        ids = set(ids)
        for record in parse_fastq(f_in):
            record_id=record[0].split(" ")[0]
            record_id=re.sub('\/1', '', record_id)
            record_id=re.sub('\/2', '', record_id)
                if record_id not in ids:
                    write_fastq(record, f_out)
                else:
                    log.write(f"{record[0]} filtered\n")		


def remove_pair_id(id: str, log: TextIO) -> str:
    """Remove the 1 or 2 from a paired read ID

    id: id string
    """
    id = id.strip()
    if id[-2:] == "/1" or id[-2:] == "/2":
        return id[:-2]

    # Assuming it's the newer id variant where komplexity removes the second half (containing pair number)
    return id
