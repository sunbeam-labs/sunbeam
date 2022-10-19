import pandas
from sunbeamlib import reports

quality_list = [reports.parse_fastqc_quality(file) for file in snakemake.input.files]
quality_table = pandas.concat(quality_list, axis=1).transpose()
quality_table.to_csv(snakemake.output[0], sep="\t", index_label="Samples")