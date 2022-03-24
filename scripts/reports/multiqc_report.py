import os

report_name = snakemake.output[0].split('/')[-1]  # Get unique name from targets.rules file
print("THIS IS THE REPORT NAME: " + report_name)
print("OUTPUT DIR: " + snakemake.params.outdir)
os.system("multiqc -i \"{a}\" -n {b} -o {c} {c}".format(
            a=snakemake.params.title, b=report_name, c=snakemake.params.outdir))