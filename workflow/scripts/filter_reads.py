import os
import subprocess
from collections import OrderedDict
from io import TextIOWrapper


def count_host_reads(fp: str, hostdict: dict, net_hostlist: set):
    hostname = os.path.basename(os.path.dirname(fp))
    hostcts = int(subprocess.getoutput("cat {} | wc -l".format(fp)).strip())
    hostdict[hostname] = hostcts

    with open(fp) as f:
        for l in f.readlines():
            net_hostlist.add(l)  # Only adds unique ids


def calculate_counts(fp: str, net_hostlist: set) -> tuple:
    original = int(str(subprocess.getoutput("zcat {} | wc -l".format(fp))).strip()) // 4
    host = len(net_hostlist)
    nonhost = int(original - host)

    return host, nonhost


def write_log(f: TextIOWrapper, hostdict: dict, host: int, nonhost: int):
    f.write("{}\n".format("\t".join(list(hostdict.keys()) + ["host", "nonhost"])))
    f.write(
        "{}\n".format("\t".join(map(str, list(hostdict.values()) + [host, nonhost])))
    )


hostdict = OrderedDict()
net_hostlist = set()
for hostid in snakemake.input.hostids:
    count_host_reads(hostid, hostdict, net_hostlist)

host, nonhost = calculate_counts(snakemake.input.reads, net_hostlist)

os.system(
    """
gzip -dc {a} | \
rbt fastq-filter {b} | \
gzip > {c}
""".format(
        a=snakemake.input.reads, b=snakemake.input.hostreads, c=snakemake.output.reads
    )
)

with open(snakemake.output.log, "w") as log:
    write_log(log, hostdict, host, nonhost)
