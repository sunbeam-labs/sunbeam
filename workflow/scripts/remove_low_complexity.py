from sunbeamlib.qc import filter_ids

ids = []
with open(snakemake.input.ids) as f:
    ids = [id.strip() for id in f.readlines()]

filter_ids(snakemake.input.reads, snakemake.output[0], ids)