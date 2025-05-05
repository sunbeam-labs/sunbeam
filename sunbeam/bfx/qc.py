"""
Supporting functions for QC rules.
"""

import gzip
import sys
from pathlib import Path
from sunbeam.bfx.parse import parse_fastq, write_fastq
from typing import Set, TextIO


def filter_ids(fp_in: Path, fp_out: Path, ids: Set[str], log: TextIO) -> None:
    """
    Filter FASTQ records based on a set of IDs to remove.

    Args:
        fp_in (Path): Path to the input FASTQ file.
        fp_out (Path): Path to the output FASTQ file.
        ids (Set[str]): Set of IDs to filter.
        log (TextIO): TextIO object to write log messages.

    Returns:
        None: This function does not return anything.

    Raises:
        AssertionError: If the number of removed IDs does not match the expected count.

    """
    with gzip.open(fp_in, "rt") as f_in, gzip.open(fp_out, "wt") as f_out:
        num_ids = len(ids)
        counter = 0
        counter_kept = 0
        for record in parse_fastq(f_in):
            counter += 1
            record = (remove_pair_id(record[0], log), record[1], record[2], record[3])
            if record[0] in ids:
                ids.remove(record[0])
            else:
                counter_kept += 1
                write_fastq(record, f_out)

        if counter - counter_kept != num_ids:
            log.write(
                f"ERROR: Mismatch (Removed: {counter - counter_kept}, Supposed to remove: {num_ids})\n"
            )
            log.write(f"IDs not found: {ids}\n")
            assert (
                False
            ), f"ERROR: Mismatch (Removed: {counter - counter_kept}, Supposed to remove: {num_ids})"

        if len(ids) > 0:
            log.write(f"WARNING: {len(ids)} ids not found in FASTQ\n")
            log.write(f"IDs not found: {ids}\n")
        else:
            log.write("IDs list empty, finished filtering\n")


def remove_pair_id(id: str, log: TextIO = sys.stdout) -> str:
    """
    Removes the pair identifier from the given ID.

    Args:
        id (str): The ID to remove the pair identifier from.
        log (TextIO): The log file to write any messages to.

    Returns:
        str: The ID with the pair identifier removed.
    """
    id = id.strip()
    if id[-2:] == "/1" or id[-2:] == "/2":
        return id[:-2]

    if " " in id:
        return id.split(" ")[0]

    return id
