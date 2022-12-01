import csv
from Bio import SearchIO
from xml.etree.ElementTree import ParseError

with open(snakemake.output[0], "w") as out:
    writer = csv.DictWriter(out, fieldnames=["sample", "query", "hit"], delimiter="\t")
    writer.writeheader()
    list(writer.writerow(result) for result in blast_summary(snakemake.input))
