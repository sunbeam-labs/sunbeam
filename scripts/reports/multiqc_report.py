import os
import shutil

report_name = snakemake.output[0].split('/')[-1]  # Get unique name from targets.rules file

os.system("conda list")
os.system("multiqc -f -i \"{a}\" -n {b} -o {c} {c}".format(
            a=snakemake.params.title, b=report_name, c=snakemake.params.outdir))
# -f overwrites any previous reports instead of iterating the next one's name
print(os.listdir(snakemake.params.outdir))

# Different versions of multiqc seem to be inconsistent with naming conventions
if report_name in os.listdir(snakemake.params.outdir):
    None
else:
    for fp in os.listdir(snakemake.params.outdir):
        if 'multiqc' in fp:
            shutil.copyfile(os.path.join(snakemake.params.outdir, fp), os.path.join(snakemake.params.outdir, report_name))
