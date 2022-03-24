import os

fwd_adapters = snakemake.config['qc']['fwd_adapters']
rev_adapters = snakemake.config['qc']['rev_adapters']
if fwd_adapters or rev_adapters:
    overlap = float('inf')
    if fwd_adapters:
        overlap = min(min(len(a) for a in fwd_adapters), overlap)
        fwd_adapter_str = "-b " + " -b ".join(snakemake.config['qc']['fwd_adapters'])
    if rev_adapters:
        overlap = min(min(len(a) for a in rev_adapters), overlap)
        rev_adapter_str = "-B " + " -B ".join(snakemake.config['qc']['rev_adapters'])
    os.system("""
    cutadapt --discard-trimmed -O {a} \
    --cores {b} \
    {c} {d} \
    -o {e} -p {f} \
    {g} {h} \
    > {j} 2>&1
    gzip {e}
    gzip {f}
    """.format(a=overlap, b=snakemake.threads, c=fwd_adapter_str, d=rev_adapter_str, e=snakemake.params.r1, f=snakemake.params.r2,
            g=snakemake.input.r1, h=snakemake.input.r2, j=snakemake.log))
else:
    os.system("""
    ln -s {a} {b} && ln -s {c} {d}
    """.format(a=snakemake.input.r1, b=snakemake.output.gr1, c=snakemake.input.r2, d=snakemake.output.gr2))