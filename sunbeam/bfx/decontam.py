import gzip
import os
import shutil
import subprocess as sp
from collections import OrderedDict
from io import TextIOWrapper
from pathlib import Path
import sys
from sunbeam.bfx.parse import parse_fastq, write_fastq


from typing import Dict, Iterator, Tuple, Union
from sunbeam.bfx.parse import parse_sam


def count_host_reads(fp: str, hostdict: dict) -> set:
    hostname = os.path.basename(os.path.dirname(fp))
    hostcts = int(sp.getoutput("cat {} | wc -l".format(fp)).strip())
    hostdict[hostname] = hostcts
    with open(fp) as f:
        return set(l.strip() for l in f.readlines())


def calculate_counts(fp: str, len_hostlist: int) -> tuple:
    with gzip.open(fp, "rt") as f:
        input_lines = sum(1 for _ in f)
    original = input_lines // 4
    nonhost = int(original - len_hostlist)
    return len_hostlist, nonhost


def write_log(f: TextIOWrapper, hostdict: OrderedDict, host: int, nonhost: int):
    f.write("{}\n".format("\t".join(list(hostdict.keys()) + ["host", "nonhost"])))
    f.write(
        "{}\n".format("\t".join(map(str, list(hostdict.values()) + [host, nonhost])))
    )


def filter_host_reads(
        input_hostids_fps, input_hostreads_fp, input_reads_fp, output_reads_fp,
        output_log_fp, log_fp):
    with open(log_fp, "w") as l:
        l.write("In filter_host_reads\n")
        hostdict = OrderedDict()
        done = False
        net_hostlist = set()
        for hostid in sorted(input_hostids_fps):
            net_hostlist.update(count_host_reads(hostid, hostdict))

        host, nonhost = calculate_counts(input_reads_fp, len(net_hostlist))

        # Check for empty host reads file
        with open(input_hostreads_fp) as f:
            # TODO: Remove aggregate_reads rule and just handle the host ids files here
            if not f.readline():
                s = f"WARNING: {input_hostreads_fp} is empty, skipping...\n"
                l.write(s)
                sys.stderr.write(s)
                shutil.copyfile(input_reads_fp, output_reads_fp)
                done = True

        # Perform filtering if host reads file is not empty
        if not done:
            with (
                    gzip.open(input_reads_fp, "rt") as f_in,
                    gzip.open(output_reads_fp, "wt") as f_out,
            ):
                for header_str, seq_str, plus_str, quality_str in parse_fastq(f_in):
                    parsed_header = (
                        header_str.split(" ")[0].replace("/1", "").replace("/2", "")
                    )
                    if not parsed_header in net_hostlist:
                        write_fastq([header_str, seq_str, plus_str, quality_str], f_out)

            # Check that the output file is about the right size given the number of ids removed
            if (
                Path(input_reads_fp).stat().st_size
                == Path(output_reads_fp).stat().st_size
                and Path(input_hostreads_fp).stat().st_size != 0
            ):
                s = f"ERROR: {input_hostreads_fp} is not empty but {input_reads_fp} and {output_reads_fp} are the same size. Something went wrong in the filtering."
                l.write(s)
                sys.stderr.write(s)
                sys.exit(1)

        with open(output_log_fp, "w") as log:
            write_log(log, hostdict, host, nonhost)

        sys.stderr.write("filter_reads script finished\n")


def get_mapped_reads(fp: str, min_pct_id: float, min_len_frac: float) -> Iterator[str]:
    """
    Takes a SAM file and returns an iterator of read names that are mapped
    """
    with open(fp, "r") as sam_file:
        for read in parse_sam(sam_file):
            if (
                (not read["FLAG"] & 0x4)  # not unmapped
                and (_get_frac(read) > min_len_frac)
                and (_get_pct_identity(read) > min_pct_id)
            ):
                yield read["QNAME"]


def _get_pct_identity(
    read: Dict[str, Union[int, float, str, Tuple[int, str]]],
) -> float:
    edit_dist = read.get("NM", 0)
    pct_mm = float(edit_dist) / len(read["SEQ"])
    return 1 - pct_mm


def _get_frac(read: Dict[str, Union[int, float, str, Tuple[int, str]]]) -> float:
    clip = 0
    for pair in read["CIGAR"]:
        if pair[1] == "S" or pair[1] == "H":
            clip += pair[0]
    frac = float(len(read["SEQ"])) / (len(read["SEQ"]) + clip)
    return frac
