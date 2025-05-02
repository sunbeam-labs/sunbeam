import gzip
from sunbeam.bfx import filter_ids, remove_pair_id


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
