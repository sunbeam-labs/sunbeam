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
        rev_adapter_str = "-g " + " -g ".join(snakemake.config['qc']['rev_adapters'])
    os.system("""
    cutadapt --discard-trimmed -O {a} \
    --cores {b} \
    {c} {d} \
    -o {e} \
    {f} \
    > {g} 2>&1
    gzip {f} || echo "Already zipped"
    ln -s {f} {h}
    """.format(a=overlap, b=snakemake.threads, c=fwd_adapter_str, d=rev_adapter_str, e=snakemake.params.tmp, f=snakemake.input,
            g=snakemake.log, h=snakemake.output))
else:
    os.system("""
    ln -s {a} {b}
    """.format(a=snakemake.input, b=snakemake.output))