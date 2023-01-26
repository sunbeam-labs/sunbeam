import subprocess as sp
import sys
from collections import OrderedDict
from filter_reads_f import count_host_reads, calculate_counts, write_log

with open(snakemake.log[0], "w") as l:
    sys.stderr = sys.stdout = l
    hostdict = OrderedDict()
    net_hostlist = set()
    for hostid in snakemake.input.hostids:
        count_host_reads(hostid, hostdict, net_hostlist)

    host, nonhost = calculate_counts(snakemake.input.reads, net_hostlist)

    sp.check_output(
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
    )

    with open(snakemake.output.log, "w") as log:
        write_log(log, hostdict, host, nonhost)
