import sys
from sunbeamlib import qc
from pathlib import Path

with open(snakemake.log[0], "w") as l:
    if snakemake.params.suffix:
        qc.strip_seq_id_suffix(
            snakemake.input[0], snakemake.output[0], snakemake.params.suffix
        )
    else:
        Path(snakemake.output[0]).symlink_to(snakemake.input[0])
