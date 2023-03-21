from itertools import groupby
from more_itertools import grouper


BLAST6_DEFAULTS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart", "qend", "sstart", "send",
    "evalue",
    "bitscore",
]


# Source: https://www.biostars.org/p/710/
def parse_fasta(f):
    faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        header_str = header.__next__()[1:].strip()

        # join all sequence lines to one.
        seq_str = "".join(s.strip() for s in faiter.__next__())

        yield (header_str, seq_str)


def write_fasta(record, f):
    f.write(f">{record[0]}\n")
    f.write(f"{record[1]}\n")


def parse_fastq(f):
    for g in grouper(f.readlines(), 4):
        header_str = g[0][1:].strip()
        seq_str = g[1].strip()
        plus_str = g[2].strip()
        quality_str = g[3].strip()

        yield (header_str, seq_str, plus_str, quality_str)


def write_fastq(record, f):
    for i, l in enumerate(record):
        if i == 0:
            f.write(f"@{l}\n")
        else:
            f.write(f"{l}\n")


def write_many_fastq(record_list, f):
    record_list = [
        [f"@{r[0]}\n", f"{r[1]}\n", f"{r[2]}\n", f"{r[3]}\n"] for r in record_list
    ]
    record_list = [item for sublist in record_list for item in sublist]
    f.writelines(record_list)


def parse_blast6(f, outfmt=BLAST6_DEFAULTS):
    for line in f.readlines():
        vals = line.strip().split("\t")
        if len(outfmt) == len(vals):
            yield dict(zip(outfmt, vals))
