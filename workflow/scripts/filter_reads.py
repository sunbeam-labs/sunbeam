import gzip
import os
import shutil
import subprocess as sp
import sys
from collections import OrderedDict
from io import TextIOWrapper
from sunbeamlib.parse import parse_fastq, write_fastq


def count_host_reads(fp: str, hostdict: dict, net_hostlist: set):
    hostname = os.path.basename(os.path.dirname(fp))
    hostcts = int(sp.getoutput("cat {} | wc -l".format(fp)).strip())
    hostdict[hostname] = hostcts

    with open(fp) as f:
        for l in f.readlines():
            net_hostlist.add(l)  # Only adds unique ids


def calculate_counts(fp: str, net_hostlist: set) -> tuple:
    original = int(str(sp.getoutput("zcat {} | wc -l".format(fp))).strip()) // 4
    host = len(net_hostlist)
    nonhost = int(original - host)

    return host, nonhost


def write_log(f: TextIOWrapper, hostdict: dict, host: int, nonhost: int):
    f.write("{}\n".format("\t".join(list(hostdict.keys()) + ["host", "nonhost"])))
    f.write(
        "{}\n".format("\t".join(map(str, list(hostdict.values()) + [host, nonhost])))
    )


with open(snakemake.log[0], "w") as l:
    hostdict = OrderedDict()
    done = False
    net_hostlist = set()
    for hostid in snakemake.input.hostids:
        count_host_reads(hostid, hostdict, net_hostlist)

    host, nonhost = calculate_counts(snakemake.input.reads, net_hostlist)

    with open(snakemake.input.hostreads) as f:
        if not f.readlines():
            s = f"WARNING: {snakemake.input.hostreads} is empty, skipping...\n"
            l.write(s)
            sys.stderr.write(s)
            shutil.copyfile(snakemake.input.reads, snakemake.output.reads)
            done = True

    if not done:
        with gzip.open(snakemake.input.reads, "rt") as f_in, open(snakemake.input.hostreads) as f_ids, gzip.open(snakemake.output.reads, "wt") as f_out:
            ids = [x.strip() for x in f_ids.readlines()]
            for header_str, seq_str, plus_str, quality_str in parse_fastq(f_in):
                if any([id in header_str for id in ids]):
                    write_fastq([header_str, seq_str, plus_str, quality_str], f_out)

    with open(snakemake.output.log, "w") as log:
        write_log(log, hostdict, host, nonhost)
    
    sys.stderr.write("filter_reads script finished\n")
