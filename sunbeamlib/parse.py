from itertools import groupby
from more_itertools import grouper

# Source: https://www.biostars.org/p/710/
def parse_fasta(f):
    faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        header_str = header.__next__()[1:].strip()

        # join all sequence lines to one.
        seq_str = "".join(s.strip() for s in faiter.__next__())

        yield (header_str, seq_str)


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