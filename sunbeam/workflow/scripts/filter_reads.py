import gzip
import os
import shutil
import subprocess as sp
import sys
from collections import OrderedDict
from io import TextIOWrapper
from pathlib import Path
from sunbeam.bfx.parse import parse_fastq, write_fastq


def count_host_reads(fp: str, hostdict: dict) -> set:
    hostname = os.path.basename(os.path.dirname(fp))
    hostcts = int(sp.getoutput("cat {} | wc -l".format(fp)).strip())
    hostdict[hostname] = hostcts

    with open(fp) as f:
        return set(l.strip() for l in f.readlines())


def calculate_counts(fp: str, len_hostlist: int) -> tuple:
    original = int(str(sp.getoutput("zcat {} | wc -l".format(fp))).strip()) // 4
    nonhost = int(original - len_hostlist)

    return len_hostlist, nonhost


def write_log(f: TextIOWrapper, hostdict: OrderedDict, host: int, nonhost: int):
    f.write("{}\n".format("\t".join(list(hostdict.keys()) + ["host", "nonhost"])))
    f.write(
        "{}\n".format("\t".join(map(str, list(hostdict.values()) + [host, nonhost])))
    )


with open(snakemake.log[0], "w") as l:
    hostdict = OrderedDict()
    done = False
    net_hostlist = set()
    for hostid in sorted(snakemake.input.hostids):
        net_hostlist.update(count_host_reads(hostid, hostdict))

    host, nonhost = calculate_counts(snakemake.input.reads, len(net_hostlist))

    # Check for empty host reads file
    with open(snakemake.input.hostreads) as f:
        # TODO: Remove aggregate_reads rule and just handle the host ids files here
        if not f.readline():
            s = f"WARNING: {snakemake.input.hostreads} is empty, skipping...\n"
            l.write(s)
            sys.stderr.write(s)
            shutil.copyfile(snakemake.input.reads, snakemake.output.reads)
            done = True

    # Perform filtering if host reads file is not empty
    if not done:
        with (
            gzip.open(snakemake.input.reads, "rt") as f_in,
            gzip.open(snakemake.output.reads, "wt") as f_out,
        ):
            for header_str, seq_str, plus_str, quality_str in parse_fastq(f_in):
                parsed_header = (
                    header_str.split(" ")[0].replace("/1", "").replace("/2", "")
                )
                if not parsed_header in net_hostlist:
                    write_fastq([header_str, seq_str, plus_str, quality_str], f_out)

        # Check that the output file is about the right size given the number of ids removed
        if (
            Path(snakemake.input.reads).stat().st_size
            == Path(snakemake.output.reads).stat().st_size
            and Path(snakemake.input.hostreads).stat().st_size != 0
        ):
            s = f"ERROR: {snakemake.input.hostreads} is not empty but {snakemake.input.reads} and {snakemake.output.reads} are the same size. Something went wrong in the filtering."
            l.write(s)
            sys.stderr.write(s)
            sys.exit(1)

    with open(snakemake.output.log, "w") as log:
        write_log(log, hostdict, host, nonhost)

    sys.stderr.write("filter_reads script finished\n")
