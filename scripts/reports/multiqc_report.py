import os

report_name = snakemake.output[0].split('/')[-1]  # Get unique name from targets.rules file

os.system("multiqc -f -i \"{a}\" -n {b} -o {c} {c}".format(
            a=snakemake.params.title, b=report_name, c=snakemake.params.outdir))
# -f overwrites any previous reports instead of iterating the next one's name

# Different versions of multiqc seem to be inconsistent with naming conventions
for fp in os.listdir(snakemake.params.outdir):
    if 'multiqc' in fp:
        os.replace(os.path.join(snakemake.params.outdir, fp), snakemake.output[0])
