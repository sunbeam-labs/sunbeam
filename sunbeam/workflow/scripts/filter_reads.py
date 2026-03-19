from sunbeam.bfx.decontam import filter_host_reads

filter_host_reads(
    input_hostids_fps=snakemake.input.hostids,
    input_hostreads_fp=snakemake.input.hostreads,
    input_reads_fp=snakemake.input.reads,
    output_reads_fp=snakemake.output.reads,
    output_log_fp=snakemake.output.log,
    log_fp=snakemake.log[0],
    )
