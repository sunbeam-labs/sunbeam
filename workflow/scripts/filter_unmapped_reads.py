import subprocess as sp
import sys
from filter_reads_f import count_host_reads, calculate_counts, write_log

with open(snakemake.log[0], "w") as l:
    net_hostlist = set()
    
    
    host_counts, nonhost_counts = calculate_counts(snakemake.input.reads, net_hostlist)