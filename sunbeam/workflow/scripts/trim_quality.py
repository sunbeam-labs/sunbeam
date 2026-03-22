from sunbeam.bfx.qc import trim_by_quality

trim_by_quality(
    input_reads_fps=snakemake.input.reads,
    output_reads_fps=snakemake.output.reads,
    output_report_fp=snakemake.output.report,
    log_fp=snakemake.log[0],
    window_width=snakemake.params.window[0],
    window_threshold=snakemake.params.window[1],
    start_threshold=snakemake.params.start_threshold,
    end_threshold=snakemake.params.end_threshold,
    min_length=snakemake.params.min_length,
    threads=snakemake.threads,
    compression=snakemake.params.compression,
)
