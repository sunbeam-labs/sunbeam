from sunbeam.bfx.reports import make_fastqc_report

make_fastqc_report(
    snakemake.input.reports,
    snakemake.output.report,
    snakemake.log[0],
    )

