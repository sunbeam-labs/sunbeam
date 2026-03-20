"""
Supporting functions for QC rules.
"""

import gzip
import sys
from pathlib import Path
import os
import subprocess as sp
from sunbeam.bfx.parse import parse_fastq, write_fastq
from typing import Set, TextIO

import traceback

def trim_by_quality(
        input_reads_fps,
        output_reads_fps,
        output_report_fp,
        log_fp,
        window_width,
        window_threshold,
        start_threshold,
        end_threshold,
        min_length,
        threads,
        compression,
        ):
    args = [
        "heyfastq",
        "trim-qual",
        "--input",
        *input_reads_fps,
        "--output",
        *output_reads_fps,
        "--report",
        output_report_fp,
        "--window-width",
        str(window_width),
        "--window-threshold",
        str(window_threshold),
        "--start-threshold",
        str(start_threshold),
        "--end-threshold",
        str(end_threshold),
        "--min-length",
        str(min_length),
        "--threads",
        str(threads),
    ]

    os.environ["HFQ_GZIP_COMPRESSION"] = str(compression)

    with open(log_fp, "w") as log:
        log.write("Initialized log and starting script...\n")
        stderr_capture = None
        try:
            res = sp.run(args, capture_output=True, text=True)
            log.write(res.stdout)
            log.write(res.stderr)
            if res.returncode != 0:
                raise RuntimeError(
                    f"heyfastq trim-qual failed with exit code {res.returncode}")
        except BaseException as e:
            log.write(f"Error during run: {e}\n")
            log.write(traceback.format_exc())
            raise
        else:
            log.write("Completed successfully.\n")


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
