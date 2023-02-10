import gzip
import os

with open(snakemake.log[0], "w") as l:
    total_count = sum(1 for line in gzip.open(snakemake.input.reads, "r")) // 4

    net_unmapped_reads = {}
    host_unmapped_counts = {}
    collisions = 0
    for i in snakemake.input.unmapped_reads:
        host_unmapped_counts[os.path.basename(os.path.dirname(i))] = total_count - sum(1 for line in open(i)) // 4
        with open(i) as f:
            id = ""
            for i, line in enumerate(f.readlines()):
                if i % 4 == 0:
                    id = line.strip()
                elif i % 4 == 1:
                    if id in net_unmapped_reads:
                        collisions += 1
                    net_unmapped_reads[id] = [line.strip()]
                else:
                    net_unmapped_reads[id].append(line.strip())
    
    total_unmapped_count = len(net_unmapped_reads)

    l.write(f"Per host unmapped count: {host_unmapped_counts}\n")
    l.write(f"Total count: {total_count}\n")
    l.write(f"Total unmapped count: {total_unmapped_count}\n")
    l.write(f"Collisions count: {collisions}\n")

    # Sanity check
    

    host = total_count - total_unmapped_count
    nonhost = total_unmapped_count
    
    with open(snakemake.output.log, "w") as f:
        f.write("{}\n".format('\t'.join(list(host_unmapped_counts.keys()) + ['host', 'nonhost'])))
        f.write(
            "{}\n".format('\t'.join(map(str, list(host_unmapped_counts.values()) + [host, nonhost])))
        )
    
    with open(snakemake.output.reads, "w") as f:
        for k, v in net_unmapped_reads.items():
            f.write(f"{k}\n")
            for val in v:
                f.write(f"{val}\n")
            
    