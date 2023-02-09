import os
import subprocess as sp
from io import TextIOWrapper


def count_host_reads(fp: str, hostdict: dict, net_hostlist: set):
    hostname = os.path.basename(os.path.dirname(fp))
    hostcts = int(sp.getoutput(f"cat {fp} | wc -l").strip())
    hostdict[hostname] = hostcts

    with open(fp) as f:
        for l in f.readlines():
            net_hostlist.add(l)  # Only adds unique ids


def calculate_counts(fp: str, net_hostlist: set) -> tuple:
    original = int(str(sp.getoutput(f"zcat {fp} | wc -l")).strip()) // 4
    host = len(net_hostlist)
    nonhost = int(original - host)

    return host, nonhost


def write_log(f: TextIOWrapper, hostdict: dict, host: int, nonhost: int):
    f.write(f"{'\t'.join(list(hostdict.keys()) + ['host', 'nonhost'])}\n")
    f.write(
        f"{'\t'.join(map(str, list(hostdict.values()) + [host, nonhost]))}\n"
    )
