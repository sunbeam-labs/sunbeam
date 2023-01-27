import subprocess as sp
import sys
from collections import OrderedDict
from filter_reads_f import count_host_reads, calculate_counts, write_log

with open(snakemake.log[0], "w") as l:
    hostdict = OrderedDict()
    net_hostlist = set()
    for hostid in snakemake.input.hostids:
        count_host_reads(hostid, hostdict, net_hostlist)

    host, nonhost = calculate_counts(snakemake.input.reads, net_hostlist)

    sys.stderr.write(f"{snakemake.input.reads}\n")
    sys.stderr.write(f"{snakemake.input.hostreads}\n")
    sys.stderr.write(f"{snakemake.output.reads}\n")
    sys.stderr.write(str([
        "gzip",
        "-dc",
        snakemake.input.reads,
        "|",
        "rbt",
        "fastq-filter",
        snakemake.input.hostreads,
        "|",
        "gzip",
        ">",
        snakemake.output.reads,
    ]))

    sp.call([
        "gzip",
        "-dc",
        snakemake.input.reads,
        "|",
        "rbt",
        "fastq-filter",
        snakemake.input.hostreads,
        "|",
        "gzip",
        ">",
        snakemake.output.reads,
    ])

    with open(snakemake.output.log, "w") as log:
        write_log(log, hostdict, host, nonhost)
