"""
Supporting functions for QC rules.
"""

import gzip
from pathlib import Path
from sunbeamlib.parse import parse_fastq, write_many_fastq
from typing import List, TextIO


def filter_ids(fp_in: Path, fp_out: Path, ids: List[str], log: TextIO) -> None:
    """Remove ids from FASTQ file.

    fp_in: path to input FASTQ
    fp_out: path to output FASTQ
    ids: list of ids to be removed
    """
    with gzip.open(fp_in, "rt") as f_in, gzip.open(fp_out, "wt") as f_out:
        records = [r for r in parse_fastq(f_in)]
        records.sort(key=lambda t: t[0])
        ids.sort()
        # Use list(records) so that it's a different object in memory and
        # you're free to remove items from the original
        for record in list(records):
            if not ids:
                log.write("IDs list empty, finished filtering\n")
                break
            if ids[0] in record[0]:
                log.write(f"{record[0]} filtered\n")
                ids.pop(0)
                records.remove(record)

        write_many_fastq(records, f_out)


def remove_pair_id(id: str, log: TextIO) -> str:
    """Remove the 1 or 2 from a paired read ID

    id: id string
    """
    id = id.strip()
    if id[-2:] == "/1" or id[-2:] == "/2":
        return id[:-2]

    # Assuming it's the newer id variant where komplexity removes the second half (containing pair number)
    return id
