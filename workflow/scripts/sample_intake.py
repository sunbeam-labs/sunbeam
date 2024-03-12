import shutil
from pathlib import Path

with open(snakemake.log[0], "w") as log:
    assert snakemake.input[0].endswith(".fastq.gz")
    log.write("Creating symlink\n")
    # Path(snakemake.output[0]).symlink_to(snakemake.input[0])
    shutil.copy(snakemake.input[0], snakemake.output[0])
