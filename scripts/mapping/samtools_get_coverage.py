from sunbeamlib import samtools

samtools.get_coverage_stats(
    snakemake.wildcards.genome,
    snakemake.input[0],
    snakemake.wildcards.sample,
    snakemake.output[0],
)
