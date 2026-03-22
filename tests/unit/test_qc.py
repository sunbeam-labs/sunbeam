import gzip
from sunbeam.bfx.qc import trim_by_quality
from sunbeam.bfx import filter_ids, remove_pair_id

input1_reads = """\
@seq
ACGTGACTGCATGACTGACTGCATACGTACGAC
+
IIIIIIIIIIIIIIIIIII++++++++++++++
"""

output1_reads = """\
@seq
ACGTGACTGCATGACTGAC
+
IIIIIIIIIIIIIIIIIII
"""

input2_reads = """\
@seq2
CAGTACGTCTCAGATCGCAGACTCAGCACGTAG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""

output2_reads = input2_reads


def test_trim_by_quality(tmp_path):
    input1_fp = tmp_path / "input1.fastq.gz"
    with gzip.open(input1_fp, "wt") as f:
        f.write(input1_reads)
    input2_fp = tmp_path / "input2.fastq.gz"
    with gzip.open(input2_fp, "wt") as f:
        f.write(input2_reads)
    output1_fp = tmp_path / "output1.fastq.gz"
    output2_fp = tmp_path / "output2.fastq.gz"
    report_fp = tmp_path / "report.txt"
    log_fp = tmp_path / "log.txt"

    trim_by_quality(
        [input1_fp, input2_fp],
        [output1_fp, output2_fp],
        report_fp,
        log_fp,
        4,
        20,
        3,
        3,
        10,
        1,
        6,
    )

    with gzip.open(output1_fp, "rt") as f:
        assert f.read() == output1_reads

    with gzip.open(output2_fp, "rt") as f:
        assert f.read() == output2_reads


def test_filter_ids(tmp_path):
    # Create a temporary input FASTQ file
    input_fastq = tmp_path / "input.fastq.gz"
    with gzip.open(input_fastq, "wt") as f:
        f.write(
            "@SEQ_ID\n"
            "GATTTGGGGTTTAAAGGGAA\n"
            "+\n"
            "IIIIIIIIIIIIIIIIII\n"
            "@ANOTHER_SEQ_ID\n"
            "GATTTGGGGTTTAAAGGGAA\n"
            "+\n"
            "IIIIIIIIIIIIIIIIII\n"
        )

    # Create a temporary output FASTQ file
    output_fastq = tmp_path / "output.fastq.gz"

    # Create a set of IDs to filter
    ids_to_filter = {"SEQ_ID"}

    # Create a log file
    log_file = tmp_path / "log.txt"
    with open(log_file, "w") as log:
        # Call the filter_ids function
        filter_ids(input_fastq, output_fastq, ids_to_filter, log)

    # Check the output FASTQ file
    with gzip.open(output_fastq, "rt") as f:
        output_content = f.read()

    assert "@SEQ_ID" not in output_content
