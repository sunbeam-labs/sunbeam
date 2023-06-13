import gzip
import os
import shutil
import subprocess as sp
import sys
from collections import OrderedDict
from io import TextIOWrapper
from sunbeamlib.parse import parse_fastq, write_many_fastq, write_fastq


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
        ids = []
        fastqs = []
        with open(snakemake.input.hostreads) as f:
            ids = [x.strip() for x in f.readlines()]
        with gzip.open(snakemake.input.reads, "rt") as f:
            for header_str, seq_str, plus_str, quality_str in parse_fastq(f):
                fastqs.append([header_str, seq_str, plus_str, quality_str])
        
        ids.sort(key=lambda x: x)
        fastqs.sort(key=lambda x: x[0])

        with open(snakemake.params.reads, "wt") as f_out:
            for record in fastqs:
                if ids[0] in record[0]:
                    ids.pop(0)
                    write_fastq(record, f_out)
                if not ids:
                    break

        with open(snakemake.params.reads) as f_in, gzip.open(
            snakemake.output.reads, "wt"
        ) as f_out:
            f_out.writelines(f_in.readlines())

    with open(snakemake.output.log, "w") as log:
        write_log(log, hostdict, host, nonhost)

    sys.stderr.write("filter_reads script finished\n")
