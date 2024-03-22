import shutil

with open(snakemake.log[0], "w") as log:
    assert snakemake.input[0].endswith(".fastq.gz")
    shutil.copy(snakemake.input[0], snakemake.output[0])
