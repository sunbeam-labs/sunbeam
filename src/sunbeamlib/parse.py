from itertools import groupby
from more_itertools import grouper
from pathlib import Path
from typing import Dict, Iterator, List, TextIO, Tuple, Union


BLAST6_DEFAULTS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


# Source: https://www.biostars.org/p/710/
def parse_fasta(f: TextIO) -> Iterator[Tuple[str, str]]:
    faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        header_str = header.__next__()[1:].strip()

        # join all sequence lines to one.
        seq_str = "".join(s.strip() for s in faiter.__next__())

        yield (header_str, seq_str)


def read_seq_ids(fasta_fp: Union[str, Path]) -> List[Tuple[str, str]]:
    """
    Return the sequence identifiers for a given fasta filepath.
    """
    with open(str(fasta_fp)) as f:
        return list(parse_fasta(f))


def write_fasta(record: Tuple[str, str], f: TextIO) -> None:
    f.write(f">{record[0]}\n")
    f.write(f"{record[1]}\n")


def parse_fastq(f: TextIO) -> Iterator[Tuple[str, str, str, str]]:
    for g in grouper(f.readlines(), 4):
        header_str = g[0][1:].strip()
        seq_str = g[1].strip()
        plus_str = g[2].strip()
        quality_str = g[3].strip()

        yield (header_str , seq_str, plus_str, quality_str)


def write_fastq(record: Tuple[str, str, str, str], f: TextIO) -> None:
    for i, l in enumerate(record):
        if i == 0:
            f.write(f"@{l}\n")
        else:
            f.write(f"{l}\n")


def write_many_fastq(record_list: List[Tuple[str, str, str, str]], f: TextIO) -> None:
    record_list = [
        [f"@{r[0]}\n", f"{r[1]}\n", f"{r[2]}\n", f"{r[3]}\n"] for r in record_list
    ]
    record_list = [item for sublist in record_list for item in sublist]
    f.writelines(record_list)


def parse_sam(
    f: TextIO,
) -> Iterator[Dict[str, Union[int, float, str, Tuple[int, str]]]]:
    for line in f:
        if line.startswith("@"):
            continue

        fields = line.strip().split("\t")
        result = {
            "QNAME": fields[0],
            "FLAG": int(fields[1]),
            "RNAME": fields[2],
            "POS": int(fields[3]),
            "MAPQ": int(fields[4]),
            "CIGAR": fields[5],
            "RNEXT": fields[6],
            "PNEXT": int(fields[7]),
            "TLEN": int(fields[8]),
            "SEQ": fields[9],
            "QUAL": fields[10],
        }

        cigar_tuples = []
        current_length = "0"

        for char in result["CIGAR"]:
            if char.isdigit():
                current_length += char
            else:
                cigar_tuples.append((int(current_length), char))
                current_length = "0"

        result["CIGAR"] = cigar_tuples

        # Parse optional fields
        optional_fields = fields[11:]
        for field in optional_fields:
            tag, data_type, value = field.split(":")
            if data_type == "i":
                value = int(value)
            elif data_type == "f":
                value = float(value)
            result[tag] = value

        yield result
