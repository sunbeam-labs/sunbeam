import os

report_name = snakemake.output[0].split('/')[-1]  # Get unique name from targets.rules file

os.system("multiqc -i \"{a}\" -n {b} -o {c} {c}".format(
            a=snakemake.params.title, b=report_name, c=snakemake.params.outdir))
if not os.path.exists(snakemake.output[0]):
    os.rename(os.path.join(snakemake.params.outdir, report_name),
            snakemake.output[0])