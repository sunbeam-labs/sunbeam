import pandas
import sys
from sunbeamlib import reports

with open(snakemake.log[0], "w") as l:
    sys.stderr = sys.stdout = l
    quality_list = [
        reports.parse_fastqc_quality(file) for file in snakemake.input.files
    ]
    quality_list = [qr for qr in quality_list if qr is not None]
    quality_table = pandas.concat(quality_list, axis=1).transpose()
    quality_table.to_csv(snakemake.output[0], sep="\t", index_label="Samples")
