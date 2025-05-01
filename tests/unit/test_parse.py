import gzip
from sunbeam.bfx.parse import (
    parse_fasta,
    parse_fastq,
    parse_sam,
    write_fasta,
    write_fastq,
)


def test_parse_fasta(tmp_path):
    fasta_fp = tmp_path / "test.fasta"
    with open(fasta_fp, "w") as f:
        f.write(">seq1\n")
        f.write("ATCG\n")
        f.write(">seq2\n")
        f.write("GCTA\n")

    expected = [("seq1", "ATCG"), ("seq2", "GCTA")]
    with open(fasta_fp) as f:
        assert list(parse_fasta(f)) == expected

    # Test with gzip
    with open(fasta_fp, "rb") as f:
        with gzip.open(fasta_fp.with_suffix(".gz"), "wb") as gz_f:
            gz_f.writelines(f)
    with gzip.open(fasta_fp.with_suffix(".gz"), "rt") as f:
        assert list(parse_fasta(f)) == expected


def test_write_fasta(tmp_path):
    fasta_fp = tmp_path / "test.fasta"
    record = ("seq1", "ATCG")
    with open(fasta_fp, "w") as f:
        write_fasta(record, f)

    with open(fasta_fp) as f:
        assert f.read() == ">seq1\nATCG\n"


def test_parse_fastq(tmp_path):
    fastq_fp = tmp_path / "test.fastq"
    with open(fastq_fp, "w") as f:
        f.write("@seq1\n")
        f.write("ATCG\n")
        f.write("+\n")
        f.write("IIII\n")
        f.write("@seq2\n")
        f.write("GCTA\n")
        f.write("+\n")
        f.write("IIII\n")

    expected = [
        ("seq1", "ATCG", "+", "IIII"),
        ("seq2", "GCTA", "+", "IIII"),
    ]
    with open(fastq_fp) as f:
        assert list(parse_fastq(f)) == expected

    # Test with gzip
    with open(fastq_fp, "rb") as f:
        with gzip.open(fastq_fp.with_suffix(".gz"), "wb") as gz_f:
            gz_f.writelines(f)
    with gzip.open(fastq_fp.with_suffix(".gz"), "rt") as f:
        assert list(parse_fastq(f)) == expected


def test_write_fastq(tmp_path):
    fastq_fp = tmp_path / "test.fastq"
    record = ("seq1", "ATCG", "+", "IIII")
    with open(fastq_fp, "w") as f:
        write_fastq(record, f)

    with open(fastq_fp) as f:
        assert f.read() == "@seq1\nATCG\n+\nIIII\n"


def test_parse_sam(tmp_path):
    sam_fp = tmp_path / "test.sam"
    with open(sam_fp, "w") as f:
        f.write("@HD\tVN:1.0\n")
        f.write("@SQ\tSN:seq1\tLN:4\n")
        f.write("seq1\t0\tseq1\t1\t255\t4M\t*\t0\t0\tATCG\tIIII\n")

    expected = [
        {
            "QNAME": "seq1",
            "FLAG": 0,
            "RNAME": "seq1",
            "POS": 1,
            "MAPQ": 255,
            "CIGAR": [(4, "M")],
            "RNEXT": "*",
            "PNEXT": 0,
            "TLEN": 0,
            "SEQ": "ATCG",
            "QUAL": "IIII",
        }
    ]
    with open(sam_fp) as f:
        assert list(parse_sam(f)) == expected