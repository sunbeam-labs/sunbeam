import os

filename = os.path.basename(snakemake.input[0])
os.system(
"""
if [ -s {a} ]
then
    vsearch --sortbylength {a} \
    --minseqlength {b} \
    --maxseqlength -1 \
    --notrunclabels \
    --output {a}.{b}f &> {c} && \
    cp {a}.{b}f {d}
else
    cp {a} {d} &> {c}
fi
""".format(a=snakemake.input, b=snakemake.params.len, c=snakemake.log, d=snakemake.output))