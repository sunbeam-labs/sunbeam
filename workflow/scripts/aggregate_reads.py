import os
import subprocess as sp

if len(snakemake.input) == 0:
    os.system("touch {0}".format(snakemake.output))
    # sp.call(
    #    [
    #        "touch",
    #        snakemake.output,
    #    ]
    # )
else:
    os.system("cat {0} > {1}".format(snakemake.input, snakemake.output))
    # sp.call(
    #    [
    #        "cat",
    #        snakemake.input,
    #        ">",
    #        snakemake.output,
    #    ]
    # )
