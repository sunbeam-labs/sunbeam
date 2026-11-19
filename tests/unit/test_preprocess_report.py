import json
import math
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


# Shaped like the JSON written by `heyfastq trim-qual --report` (heyfastq
# 0.2.3). For paired-end input heyfastq counts read pairs and takes bases from
# R1 only, so 2000 pairs of 70 bp reads (4000 reads in fastp's terms) are
# 140000 bases.
QUALITY_EXAMPLE = {
    "version": "0.2.3",
    "chunk_size": 1000,
    "window_width": 4,
    "window_threshold": 15.0,
    "start_threshold": 3.0,
    "end_threshold": 3.0,
    "min_length": 36,
    "trim_avg": {
        "input_reads": 2000,
        "input_bases": 140000,
        "output_reads": 2000,
        "output_bases": 136100,
    },
    "trim_ends": {
        "input_reads": 2000,
        "input_bases": 136100,
        "output_reads": 2000,
        "output_bases": 136000,
    },
    "filter_length": {
        "input_reads": 2000,
        "input_bases": 136000,
        "output_reads": 1998,
        "output_bases": 135864,
    },
}

# Shaped like the JSON written by `heyfastq filter-kscore --report`.
COMPLEXITY_EXAMPLE = {
    "version": "0.2.3",
    "chunk_size": 1000,
    "kmer_size": 4,
    "min_kscore": 0.55,
    "filter_kscore": {
        "input_reads": 1998,
        "input_bases": 135864,
        "output_reads": 1995,
        "output_bases": 135660,
    },
}

# Plain-text report format that older versions of the parser expected.
LEGACY_QUALITY_EXAMPLE = """Input reads: total=4000, average_length=70.00\nOutput reads: total=4000, average_length=70.00\n"""

LEGACY_COMPLEXITY_EXAMPLE = """Input reads: total=4000, average_length=70.00\nOutput reads: total=3990, average_length=70.00\n"""

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
    report_path = write_json(tmp_path / "02_quality_sample1.json", QUALITY_EXAMPLE)
    parsed = parse_quality_report(report_path)

    # Input is before quality trimming (trim_avg), output is after the
    # minimum-length filter (filter_length)
    assert parsed == {
        "quality_input_reads": 2000,
        "quality_output_reads": 1998,
        "quality_input_average_length": pytest.approx(70),
        "quality_output_average_length": pytest.approx(68),
    }


def test_parse_complexity_report(tmp_path):
    report_path = write_json(
        tmp_path / "03_complexity_sample1.json", COMPLEXITY_EXAMPLE
    )
    parsed = parse_complexity_report(report_path)

    assert parsed == {
        "complexity_input_reads": 1998,
        "complexity_output_reads": 1995,
        "complexity_input_average_length": pytest.approx(68),
        "complexity_output_average_length": pytest.approx(68),
    }


def test_parse_quality_report_sums_multiple_reports(tmp_path):
    second = json.loads(json.dumps(QUALITY_EXAMPLE))
    second["trim_avg"].update(input_reads=1000, input_bases=50000)
    second["filter_length"].update(output_reads=500, output_bases=20000)
    report_paths = [
        write_json(tmp_path / "a.json", QUALITY_EXAMPLE),
        write_json(tmp_path / "b.json", second),
    ]
    parsed = parse_quality_report(report_paths)

    assert parsed["quality_input_reads"] == 3000
    assert parsed["quality_output_reads"] == 2498
    assert parsed["quality_input_average_length"] == pytest.approx(190000 / 3000)
    assert parsed["quality_output_average_length"] == pytest.approx(155864 / 2498)


def test_parse_complexity_report_no_reads(tmp_path):
    empty = json.loads(json.dumps(COMPLEXITY_EXAMPLE))
    empty["filter_kscore"] = dict.fromkeys(empty["filter_kscore"], 0)
    report_path = write_json(tmp_path / "complexity.json", empty)
    parsed = parse_complexity_report(report_path)

    assert parsed["complexity_input_reads"] == 0
    assert parsed["complexity_output_reads"] == 0
    assert math.isnan(parsed["complexity_input_average_length"])
    assert math.isnan(parsed["complexity_output_average_length"])


def test_parse_quality_report_rejects_unrecognised_json(tmp_path):
    # e.g. a complexity report passed where a quality report was expected
    report_path = write_json(tmp_path / "complexity.json", COMPLEXITY_EXAMPLE)

    with pytest.raises(ValueError, match="trim_avg"):
        parse_quality_report(report_path)


def test_parse_quality_report_legacy_text(tmp_path):
    report_path = write_file(tmp_path / "quality.txt", LEGACY_QUALITY_EXAMPLE)
    parsed = parse_quality_report(report_path)

    assert parsed["quality_input_reads"] == 4000
    assert parsed["quality_output_reads"] == 4000
    assert parsed["quality_input_average_length"] == pytest.approx(70)
    assert parsed["quality_output_average_length"] == pytest.approx(70)


def test_parse_complexity_report_legacy_text(tmp_path):
    report_path = write_file(tmp_path / "complexity.txt", LEGACY_COMPLEXITY_EXAMPLE)
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
sample1	4000	280000	4000	280000	4000	0.0	70.0	70.0	0.48755699999999996	0.48755699999999996	2000	1998	70.0	68.0	1998	1995	68.0	68.0	0	0	0	0	1995
"""


def test_make_preprocess_report(tmp_path):
    adapter_fp = write_json(tmp_path / "adapter.json", ADAPTER_EXAMPLE)
    trim_fp = write_json(tmp_path / "02_quality_sample1.json", QUALITY_EXAMPLE)
    complexity_fp = write_json(
        tmp_path / "03_complexity_sample1.json", COMPLEXITY_EXAMPLE
    )
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
