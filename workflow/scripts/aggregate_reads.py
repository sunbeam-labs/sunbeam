import os
import sys

if not snakemake.params.host_fp:
    with open(snakemake.log[0], "w") as f_log:
        f_log.write("NO HOST_FP SPECIFIED")
        sys.exit(1)

if len(snakemake.input) == 0:
    os.system("touch {0}".format(snakemake.output))
else:
    os.system("cat {0} > {1}".format(snakemake.input, snakemake.output))
