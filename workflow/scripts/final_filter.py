import shutil
import sys
from final_filter_f import parse_fasta, filter_seqs, write_fasta

with open(snakemake.log[0], "w") as l:
    with open(snakemake.input[0]) as f:
        with open(f"{snakemake.input[0]}.{snakemake.params.len}f", "w") as g:
            write_fasta(g, list(filter_seqs(parse_fasta(f), snakemake.params.len)))

    shutil.copyfile(
        f"{snakemake.input[0]}.{snakemake.params.len}f", snakemake.output[0]
    )
