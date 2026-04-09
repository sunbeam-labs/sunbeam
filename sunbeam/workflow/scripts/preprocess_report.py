from sunbeam.bfx.reports import make_preprocess_report

make_preprocess_report(
    adapter_fps=snakemake.input.adapter,
    trim_fps=snakemake.input.trim,
    complexity_fps=snakemake.input.complexity,
    decontam_fps=snakemake.input.decontam,
    output_fp=snakemake.output.report,
    samples=snakemake.params.samples,
    log_fp=snakemake.log[0],
)
