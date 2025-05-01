from sunbeam.bfx.qc import filter_ids, remove_pair_id

with open(snakemake.log[0], "w") as log:
    with open(snakemake.input.ids) as f:
        ids = set(remove_pair_id(id, log) for id in f.readlines())
    log.write(f"Num Komplexity IDs to be filtered: {len(ids)}\n")

    filter_ids(snakemake.input.reads, snakemake.output[0], ids, log)
