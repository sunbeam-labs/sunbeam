import gzip

from sunbeam.bfx import get_mapped_reads
from sunbeam.bfx.decontam import (
    filter_host_reads,
)


filter_host_input_fastq = """\
@NZ_CP069563.1_58552_59002_1:0:0_1:0:0_0/1
ACTGTTAAAACAAGTCTTTTGTCTTGAAAGAACAAGTCTTTTAGAAGCAAAAGCTTTGTCCTTTCCATCT
+
2222222222222222222222222222222222222222222222222222222222222222222222
@NZ_CP069563.1_898_1404_4:0:0_1:0:0_1/1
TGCGCCTGTTTCCAATGTTGTGTGCCTTGTCGTTTCTCTCGAAAACGATATATACGTGAGGAGTTATAGT
+
2222222222222222222222222222222222222222222222222222222222222222222222
@NZ_CP069563.1_23599_24082_1:0:0_1:0:0_2/1
GTCATCCGGTGGCCCTGGTTTGTAGCCTCGGTCATCAACTGCCTGGCAGGAGCATGGCTACACCTCCGGC
+
2222222222222222222222222222222222222222222222222222222222222222222222
@NZ_CP069563.1_16547_17065_1:0:0_1:0:0_3/1
TTCTCTCGAAAGCGGTTTGGGACGCTCACGCACCAGGGCAAACAGCTGTCCGGACAACTCTTTAAGCAAT
+
2222222222222222222222222222222222222222222222222222222222222222222222
"""

filter_host_input_hostids = """\
NZ_CP069563.1_58552_59002_1:0:0_1:0:0_0
NZ_CP069563.1_898_1404_4:0:0_1:0:0_1
"""

expected_filter_host_output = """\
@NZ_CP069563.1_23599_24082_1:0:0_1:0:0_2/1
GTCATCCGGTGGCCCTGGTTTGTAGCCTCGGTCATCAACTGCCTGGCAGGAGCATGGCTACACCTCCGGC
+
2222222222222222222222222222222222222222222222222222222222222222222222
@NZ_CP069563.1_16547_17065_1:0:0_1:0:0_3/1
TTCTCTCGAAAGCGGTTTGGGACGCTCACGCACCAGGGCAAACAGCTGTCCGGACAACTCTTTAAGCAAT
+
2222222222222222222222222222222222222222222222222222222222222222222222
"""

expected_filter_host_log = """\
fakehost	host	nonhost
2	2	2
"""


def test_filter_host_reads(tmp_path):
    input_fastq_fp = tmp_path / "input.fastq.gz"
    with gzip.open(input_fastq_fp, "wt") as f:
        f.write(filter_host_input_fastq)
    # column in log file derived from dirname of ID list
    fakehost_dir = tmp_path / "fakehost"
    fakehost_dir.mkdir()
    input_hostids_fp = fakehost_dir / "input.ids"
    with open(input_hostids_fp, "wt") as f:
        f.write(filter_host_input_hostids)
    input_hostreads_fp = input_hostids_fp

    output_fastq_fp = tmp_path / "output.fastq.gz"
    output_log_fp = tmp_path / "output.log"
    log_fp = tmp_path / "log.txt"
    filter_host_reads(
        [input_hostids_fp],
        input_hostreads_fp,
        input_fastq_fp,
        output_fastq_fp,
        output_log_fp,
        log_fp)

    with gzip.open(output_fastq_fp, "rt") as f:
        assert f.read() == expected_filter_host_output

    with open(output_log_fp) as f:
        assert f.read() == expected_filter_host_log


def test_get_mapped_reads(tmp_path):
    # Create a temporary SAM file
    sam_file = tmp_path / "test.sam"
    with open(sam_file, "w") as f:
        f.write(
            "@HD\tVN:1.0\tSO:unsorted\n"
            "@SQ\tSN:seq1\tLN:1000\n"
            "@SQ\tSN:seq2\tLN:2000\n"
            # Basic mapped read (perfect match)
            "read1\t0\tseq1\t1\t255\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\n"
            # Unmapped read
            "read2\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n"
            # Reverse strand mapped read (perfect match)
            "read3\t16\tseq1\t5\t255\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\n"
            # Low percent identity (3 mismatches out of 10 = 70% identity)
            "read4\t0\tseq1\t10\t255\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:3\n"
            # Partially aligned read with soft clipping (5/15 bases mapped = 33%)
            "read5\t0\tseq1\t20\t255\t5S5M5S\t*\t0\t0\tTTTTTACGTATTTTT\tIIIIIIIIIIIIIII\tNM:i:0\n"
            # Read with decent coverage (8/10 = 80%)
            "read6\t0\tseq1\t30\t255\t8M2S\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\n"
            # Read with insertions and deletions (should have good alignment % but low identity)
            "read7\t0\tseq1\t40\t255\t3M2I3M2D2M\t*\t0\t0\tACGTTACGTAC\tIIIIIIIIIII\tNM:i:4\n"
            # Secondary alignment (function doesn't filter these)
            "read8\t256\tseq1\t50\t255\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\n"
            # Properly paired read (perfect match)
            "read9\t3\tseq1\t60\t255\t10M\t=\t160\t100\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\n"
            # Supplementary alignment (function doesn't filter these)
            "read10\t2048\tseq1\t70\t255\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\n"
        )

    # Call the function with default thresholds (80% alignment, 80% identity)
    mapped_reads = list(get_mapped_reads(str(sam_file), 0.8, 0.8))

    # Check the output based on how get_mapped_reads actually works:
    # - It filters out unmapped reads (0x4 flag)
    # - It filters reads with low alignment fraction (calculation based on CIGAR)
    # - It filters reads with low percent identity (1 - NM/read_length)
    assert "read1" in mapped_reads  # perfect match
    assert "read2" not in mapped_reads  # unmapped
    assert "read3" in mapped_reads  # perfect match, reverse strand
    assert "read4" not in mapped_reads  # fails percent identity (70% < 80%)
    assert "read5" not in mapped_reads  # fails length fraction (33% < 80%)
    assert "read6" in mapped_reads  # passes both (80% alignment, 100% identity)
    assert "read7" not in mapped_reads  # fails percent identity due to NM tag
    assert (
        "read8" in mapped_reads
    )  # secondary alignment, but function doesn't filter these
    assert "read9" in mapped_reads  # paired read, perfect match
    assert (
        "read10" in mapped_reads
    )  # supplementary alignment, but function doesn't filter these

    # Expected count: 6 reads (read1, read3, read6, read8, read9, read10)
    assert len(mapped_reads) == 6

    # Test with lower thresholds to include more reads
    lenient_mapped_reads = list(get_mapped_reads(str(sam_file), 0.6, 0.3))
    print(lenient_mapped_reads)
    assert "read2" not in lenient_mapped_reads  # Still unmapped
    assert "read4" in lenient_mapped_reads  # Now should include low identity
    assert "read7" in lenient_mapped_reads  # Now passes identity threshold
    assert len(lenient_mapped_reads) == 9  # All but read2 and read5
