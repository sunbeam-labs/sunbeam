import shutil
import os

with open(snakemake.log[0], "w") as log:
    assert snakemake.input[0].endswith(".fastq.gz")
    os.symlink(snakemake.input[0], snakemake.output[0])
#    shutil.copy(snakemake.input[0], snakemake.output[0])
