import os

report_name = snakemake.output[0].split('/')[-1]  # Get unique name from targets.rules file

os.system("multiqc -f -i \"{a}\" -n {b} -o {c} {c}".format(
            a=snakemake.params.title, b=report_name, c=snakemake.params.outdir))
# -f overwrites any previous reports instead of iterating the next one's name
