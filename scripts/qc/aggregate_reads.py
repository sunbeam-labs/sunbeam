import os

if len(snakemake.input) == 0:
    os.system("touch {0}".format(snakemake.output))
else:
    os.system("cat {0} > {1}".format(snakemake.input, snakemake.output))