import pandas
from sunbeam.bfx.reports import summarize_qual_decontam

paired_end = snakemake.config["all"]["paired_end"]
summary_list = [
    summarize_qual_decontam(q, d, k, paired_end)
    for q, d, k in zip(
        snakemake.input.trim_files,
        snakemake.input.decontam_files,
        snakemake.input.komplexity_files,
    )
]
_reports = pandas.concat(summary_list)
_reports.to_csv(snakemake.output[0], sep="\t", index_label="Samples")
