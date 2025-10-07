import json
from pathlib import Path

import pandas
import pytest

from sunbeam.workflow.scripts.preprocess_report import (
    parse_adapter_report,
    parse_complexity_report,
    parse_decontam_report,
    parse_quality_report,
    parse_reports,
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


def test_parse_reports_integration(tmp_path):
    adapter = write_json(tmp_path / "adapter.json", ADAPTER_EXAMPLE)
    quality = write_file(tmp_path / "quality.txt", QUALITY_EXAMPLE)
    complexity = write_file(tmp_path / "complexity.txt", COMPLEXITY_EXAMPLE)
    decontam = write_file(tmp_path / "decontam.txt", DECONTAM_EXAMPLE)

    df = parse_reports(
        {
            "sample1": {
                "adapter": adapter,
                "trim": quality,
                "complexity": complexity,
                "decontam": decontam,
            }
        }
    )

    assert isinstance(df, pandas.DataFrame)
    assert list(df.index) == ["sample1"]
    assert df.loc["sample1", "adapter_total_reads_before"] == 4000
    assert df.loc["sample1", "quality_output_reads"] == 4000
    assert df.loc["sample1", "complexity_output_reads"] == 3990
    assert df.loc["sample1", "decontam_nonhost"] == 1995
