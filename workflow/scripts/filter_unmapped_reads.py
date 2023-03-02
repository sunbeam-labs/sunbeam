import gzip
import os
from sunbeamlib.parse import parse_fastq, write_fastq

with open(snakemake.log[0], "w") as l:
    total_count = sum(1 for line in gzip.open(snakemake.input.reads, "r")) // 4
    host_mapped_counts = {}
    unmapped_reads = {}

    for i in snakemake.input.unmapped_reads:
        basename = os.path.basename(os.path.dirname(i))
        host_mapped_counts[basename] = (
            total_count - sum(1 for line in open(i)) // 4
        )
        unmapped_reads[basename] = set()
        with open(i) as f:
            for record in parse_fastq(f):
                unmapped_reads[basename].add(record)
            
    final_unmapped_reads = set.intersection(*unmapped_reads.values())

    l.write(f"Per host mapped count: {host_mapped_counts}\n")
    l.write(f"Total count: {total_count}\n")

    host = sum(host_mapped_counts.values())
    nonhost = len(final_unmapped_reads)
    l.write(f"Total host mapped count: {host}\n")
    l.write(f"Total unmapped count: {nonhost}\n")
    l.write(f"Sanity check (host + nonhost = total): {host} + {nonhost} = {total_count}\n")
    assert host + nonhost == total_count

    if len(snakemake.input.unmapped_reads) == 0:
        l.write(f"No unmapped reads files, there are probably no host files\n")
        host = 0
        nonhost = total_count

    with open(snakemake.output.log, "w") as f:
        f.write(
            "{}\n".format(
                "\t".join(list(host_mapped_counts.keys()) + ["host", "nonhost"])
            )
        )
        f.write(
            "{}\n".format(
                "\t".join(
                    map(str, list(host_mapped_counts.values()) + [host, nonhost])
                )
            )
        )

    with gzip.open(snakemake.output.reads, "wt") as f_out:
        for record in final_unmapped_reads:
            write_fastq(record, f_out)
