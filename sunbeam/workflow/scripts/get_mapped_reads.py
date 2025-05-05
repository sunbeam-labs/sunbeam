from sunbeam.bfx.decontam import get_mapped_reads

with open(snakemake.output.ids, "w") as out:
    last = None
    for read_id in get_mapped_reads(
        snakemake.input[0], snakemake.params.pct_id, snakemake.params.frac
    ):
        if read_id == last:
            continue
        else:
            out.write(read_id + "\n")
            last = read_id
