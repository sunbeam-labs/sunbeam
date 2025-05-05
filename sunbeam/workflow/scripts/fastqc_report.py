import pandas
from sunbeam.bfx.reports import parse_fastqc_quality

with open(snakemake.log[0], "w") as l:
    quality_list = [parse_fastqc_quality(file) for file in snakemake.input.files]
    quality_list = [qr for qr in quality_list if qr is not None]
    quality_table = pandas.concat(quality_list, axis=1).transpose()
    quality_table.to_csv(snakemake.output[0], sep="\t", index_label="Samples")
