from sunbeamlib.qc import filter_ids

with open(snakemake.log[0], "w") as log:
    ids = []
    with open(snakemake.input.ids) as f:
        ids = [id.strip()[:-1] for id in f.readlines()]
    log.write(f"Komplexity IDs to be filtered: {str(ids)}\n")

    filter_ids(snakemake.input.reads, snakemake.output[0], ids, log)
