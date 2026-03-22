import json
from pathlib import Path
import pytest

from sunbeam.bfx.reports import (
    parse_adapter_report,
    parse_complexity_report,
    parse_decontam_report,
    parse_quality_report,
    make_preprocess_report,
)

ADAPTER_EXAMPLE = {
    "summary": {
        "fastp_version": "0.22.0",
        "sequencing": "paired end (70 cycles + 70 cycles)",
        "before_filtering": {
            "total_reads": 4000,
            "total_bases": 280000,
            "q20_bases": 0,
            "q30_bases": 0,
            "q20_rate": 0,
            "q30_rate": 0,
            "read1_mean_length": 70,
            "read2_mean_length": 70,
            "gc_content": 0.487557,
        },
        "after_filtering": {
            "total_reads": 4000,
            "total_bases": 280000,
            "q20_bases": 0,
            "q30_bases": 0,
            "q20_rate": 0,
            "q30_rate": 0,
            "read1_mean_length": 70,
            "read2_mean_length": 70,
            "gc_content": 0.487557,
        },
    },
    "filtering_result": {
        "passed_filter_reads": 4000,
        "low_quality_reads": 0,
        "too_many_N_reads": 0,
        "low_complexity_reads": 0,
        "too_short_reads": 0,
        "too_long_reads": 0,
    },
    "duplication": {"rate": 0},
}


QUALITY_EXAMPLE = """Input reads: total=4000, average_length=70.00\nOutput reads: total=4000, average_length=70.00\n"""

COMPLEXITY_EXAMPLE = """Input reads: total=4000, average_length=70.00\nOutput reads: total=3990, average_length=70.00\n"""

DECONTAM_EXAMPLE = "human\thuman_copy\tphix174\thost\tnonhost\n0\t0\t0\t0\t1995\n"


def write_file(path: Path, contents: str) -> str:
    path.write_text(contents)
    return str(path)


def write_json(path: Path, payload) -> str:
    path.write_text(json.dumps(payload))
    return str(path)


def test_parse_adapter_report(tmp_path):
    report_path = write_json(tmp_path / "adapter.json", ADAPTER_EXAMPLE)
    parsed = parse_adapter_report(report_path)

    assert parsed["adapter_total_reads_before"] == 4000
    assert parsed["adapter_total_reads_after"] == 4000
    assert parsed["adapter_total_bases_before"] == 280000
    assert parsed["adapter_total_bases_after"] == 280000
    assert parsed["adapter_passed_filter_reads"] == 4000
    assert parsed["adapter_gc_content_before"] == pytest.approx(0.487557)
    assert parsed["adapter_gc_content_after"] == pytest.approx(0.487557)
    assert parsed["adapter_read1_mean_length"] == pytest.approx(70)
    assert parsed["adapter_read2_mean_length"] == pytest.approx(70)


def test_parse_quality_report(tmp_path):
    report_path = write_file(tmp_path / "quality.txt", QUALITY_EXAMPLE)
    parsed = parse_quality_report(report_path)

    assert parsed["quality_input_reads"] == 4000
    assert parsed["quality_output_reads"] == 4000
    assert parsed["quality_input_average_length"] == pytest.approx(70)
    assert parsed["quality_output_average_length"] == pytest.approx(70)


def test_parse_complexity_report(tmp_path):
    report_path = write_file(tmp_path / "complexity.txt", COMPLEXITY_EXAMPLE)
    parsed = parse_complexity_report(report_path)

    assert parsed["complexity_input_reads"] == 4000
    assert parsed["complexity_output_reads"] == 3990


def test_parse_decontam_report(tmp_path):
    report_path = write_file(tmp_path / "decontam.txt", DECONTAM_EXAMPLE)
    parsed = parse_decontam_report(report_path)

    assert parsed == {
        "decontam_human": 0,
        "decontam_human_copy": 0,
        "decontam_phix174": 0,
        "decontam_host": 0,
        "decontam_nonhost": 1995,
    }


expected_preprocess_report = """\
sample	adapter_total_reads_before	adapter_total_bases_before	adapter_total_reads_after	adapter_total_bases_after	adapter_passed_filter_reads	adapter_duplication_rate	adapter_read1_mean_length	adapter_read2_mean_length	adapter_gc_content_before	adapter_gc_content_after	quality_input_reads	quality_output_reads	quality_input_average_length	quality_output_average_length	complexity_input_reads	complexity_output_reads	complexity_input_average_length	complexity_output_average_length	decontam_human	decontam_human_copy	decontam_phix174	decontam_host	decontam_nonhost
sample1	4000	280000	4000	280000	4000	0.0	70.0	70.0	0.48755699999999996	0.48755699999999996	4000	4000	70.0	70.0	4000	3990	70.0	70.0	0	0	0	0	1995
"""


def test_make_preprocess_report(tmp_path):
    adapter_fp = write_json(tmp_path / "adapter.json", ADAPTER_EXAMPLE)
    trim_fp = write_file(tmp_path / "quality.txt", QUALITY_EXAMPLE)
    complexity_fp = write_file(tmp_path / "complexity.txt", COMPLEXITY_EXAMPLE)
    decontam_fp = write_file(tmp_path / "decontam.txt", DECONTAM_EXAMPLE)
    output_fp = tmp_path / "output.tsv"
    samples = ["sample1"]
    log_fp = tmp_path / "log.txt"

    make_preprocess_report(
        [adapter_fp],
        [trim_fp],
        [complexity_fp],
        [decontam_fp],
        output_fp,
        ["sample1"],
        log_fp,
    )

    with open(output_fp) as f:
        assert f.read() == expected_preprocess_report
