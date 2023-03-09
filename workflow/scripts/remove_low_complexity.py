import gzip
import shutil
from sunbeamlib.qc import filter_ids, remove_pair_id

with open(snakemake.log[0], "w") as log:
    ids = []
    with open(snakemake.input.ids) as f:
        ids = [remove_pair_id(id, log) for id in f.readlines()]
    log.write(f"Komplexity IDs to be filtered: {str(ids)}\n")

    with gzip.open(snakemake.input.reads, "rb") as f_in, open(snakemake.input.reads[:-3], "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    filter_ids(snakemake.input.reads[:-3], snakemake.output.unzip, ids, log)

    with open(snakemake.output.unzip, "rb") as f_in, gzip.open(
        snakemake.output.out, "wb"
    ) as f_out:
        f_out.writelines(f_in.readlines())
