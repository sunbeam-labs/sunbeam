"""
Summarize stats for coverage data for each sample 
"""

from get_coverage_f import parse_depth, get_cov_stats, write_csv

input_fp = snakemake.input[0]
sample = snakemake.wildcards.sample
output_fp = snakemake.output[0]

with open(input_fp) as f:
    with open(output_fp, "w") as g:
        write_csv(g, get_cov_stats(parse_depth(f), sample))
