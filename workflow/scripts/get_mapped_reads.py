import sys
from sunbeamlib.decontam import get_mapped_reads

with open(snakemake.log[0], "w") as l:
    with open(snakemake.output.ids, "w") as out:
        last = None

        for read_id in get_mapped_reads(
            snakemake.input[0], snakemake.params.pct_id, snakemake.params.frac, l
        ):
            l.write(f"DEBUG: Processing read ID {read_id}")
            if read_id == last:
                continue
            else:
                l.write(f"DEBUG: Writing read ID {read_id}")
                out.write(read_id + "\n")
                last = read_id
