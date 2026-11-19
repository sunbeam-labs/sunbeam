import csv
import pytest
import subprocess as sp
from pathlib import Path
from sunbeam.scripts.config import main as Config
from sunbeam.scripts.init import main as Init
from sunbeam.scripts.run import main as Run


def assert_preprocess_summary_has_qc_counts(sunbeam_output: Path):
    """Every sample should have non-zero quality and complexity read counts,
    i.e. the heyfastq reports were actually parsed."""

    summary_fp = sunbeam_output / "qc" / "reports" / "preprocess_summary.tsv"
    with open(summary_fp, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert rows
    for row in rows:
        for step in ("quality", "complexity"):
            input_reads = int(row[f"{step}_input_reads"])
            output_reads = int(row[f"{step}_output_reads"])
            assert input_reads > 0, (row["sample"], step)
            assert 0 <= output_reads <= input_reads, (row["sample"], step)


def test_sunbeam_run(tmp_path, DATA_DIR):
    project_dir = tmp_path / "test"

    Init(
        [
            str(project_dir),
            "--data_fp",
            str(DATA_DIR / "reads"),
        ]
    )

    with pytest.raises(SystemExit) as excinfo:
        Run(
            [
                "--profile",
                str(project_dir),
                "--exclude",
                "all",
                "--show-failed-logs",
            ]
        )
    assert excinfo.value.code == 0

    sunbeam_output = project_dir / "sunbeam_output"
    assert sunbeam_output.exists()
    assert_preprocess_summary_has_qc_counts(sunbeam_output)


def test_sunbeam_run_with_single_end_reads(tmp_path, DATA_DIR):
    project_dir = tmp_path / "test"

    Init(
        [
            str(project_dir),
            "--data_fp",
            str(DATA_DIR / "single_end_reads"),
            "--single_end",
        ]
    )

    with pytest.raises(SystemExit) as excinfo:
        Run(
            [
                "--profile",
                str(project_dir),
                "--exclude",
                "all",
                "--show-failed-logs",
            ]
        )
    assert excinfo.value.code == 0

    sunbeam_output = project_dir / "sunbeam_output"
    assert sunbeam_output.exists()


def test_sunbeam_run_with_dirty_reads(tmp_path, DATA_DIR):
    project_dir = tmp_path / "test"

    Init(
        [
            str(project_dir),
            "--data_fp",
            str(DATA_DIR / "dirty_reads"),
        ]
    )

    with pytest.raises(SystemExit) as excinfo:
        Run(
            [
                "--profile",
                str(project_dir),
                "--exclude",
                "all",
            ]
        )
    assert excinfo.value.code == 0

    sunbeam_output = project_dir / "sunbeam_output"
    assert sunbeam_output.exists()
    assert_preprocess_summary_has_qc_counts(sunbeam_output)


def test_sunbeam_run_with_target_after_exclude(tmp_path, DATA_DIR, capsys):
    project_dir = tmp_path / "test"

    Init(
        [
            str(project_dir),
            "--data_fp",
            str(DATA_DIR / "reads"),
        ]
    )

    ret = sp.run(
        [
            "sunbeam",
            "run",
            "--profile",
            str(project_dir),
            "--exclude",
            "all",
            "clean_qc",
            "-n",
            "--quiet",
        ],
        check=True,
        capture_output=True,
    )

    assert ret.returncode == 0
    assert "clean_qc" in ret.stderr.decode("utf-8")
    assert "filter_reads" not in ret.stderr.decode("utf-8")
