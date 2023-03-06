import os
import pytest
import shutil
import subprocess as sp
import sys
from pathlib import Path

test_dir = Path(__file__).parent.parent.resolve()
sys.path.append(test_dir)
from config_fixture import output_dir, config


@pytest.fixture
def init(output_dir):
    output_dir = output_dir / "sunbeam_run"

    sp.check_output(
        [
            "sunbeam",
            "init",
            "--data_fp",
            f"{test_dir / 'data' / 'reads'}",
            output_dir,
        ]
    )

    config_str = f"qc: {{host_fp: {test_dir / 'data' / 'hosts'}}}"

    sp.check_output(
        [
            "sunbeam",
            "config",
            "modify",
            "-i",
            "-s",
            f"{config_str}",
            f"{output_dir / 'sunbeam_config.yml'}",
        ]
    )

    yield output_dir

    if os.environ.get("CI", False):
        try:
            shutil.copytree(output_dir, "output_sunbeam_run/")
        except FileExistsError as e:
            pass


def test_sunbeam_run_all(init):
    output_dir = init
    sunbeam_output_dir = output_dir / "sunbeam_output"

    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
        ]
    )

    with open(test_dir / "targets.txt") as f:
        for line in f.readlines():
            if not line.strip():
                continue
            target = sunbeam_output_dir / line.strip()
            if not target.exists():
                raise SystemExit(f"Target '{target}' not found")
            elif target.stat().st_size == 0:
                raise SystemExit(f"Target '{target}' is empty")
            else:
                print(f"Found target '{target}'")

    with open(sunbeam_output_dir / "qc" / "reports" / "preprocess_summary.tsv") as f:
        f.readline()  # Headers
        stats = f.readline().split("\t")
        assert int(stats[1]) == 400  # Input reads
        assert int(stats[2]) == 400  # Both kept by cutadapt
        assert (
            int(stats[6]) + int(stats[8]) + int(stats[11]) == 200
        )  # Human + phiX + komplexity
        assert int(stats[6]) == int(stats[7])  # Human = human_copy
        assert int(stats[9]) + int(stats[11]) == 200  # Host + komplexity
        assert int(stats[10]) == 200  # Nonhost


@pytest.fixture
def init_single_end(output_dir):
    output_dir = output_dir / "sunbeam_run_single_end"

    sp.check_output(
        [
            "sunbeam",
            "init",
            "--data_fp",
            f"{test_dir / 'data' / 'single_end_reads'}",
            "--single_end",
            "--format",
            "{sample}_R{rp}.fastq.gz",
            output_dir,
        ]
    )

    config_str = f"qc: {{host_fp: {test_dir / 'data' / 'hosts'}}}"

    sp.check_output(
        [
            "sunbeam",
            "config",
            "modify",
            "-i",
            "-s",
            f"{config_str}",
            f"{output_dir / 'sunbeam_config.yml'}",
        ]
    )

    yield output_dir

    if os.environ.get("CI", False):
        try:
            shutil.copytree(output_dir, "output_sunbeam_run_single_end/")
        except FileExistsError as e:
            pass


def test_sunbeam_run_all_single_end(init_single_end):
    output_dir = init_single_end
    sunbeam_output_dir = output_dir / "sunbeam_output"

    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
        ]
    )

    with open(test_dir / "targets_single_end.txt") as f:
        for line in f.readlines():
            if not line.strip():
                continue
            target = sunbeam_output_dir / line.strip()
            if not target.exists():
                raise SystemExit(f"Target '{target}' not found")
            elif target.stat().st_size == 0:
                raise SystemExit(f"Target '{target}' is empty")
            else:
                print(f"Found target '{target}'")

    with open(sunbeam_output_dir / "qc" / "reports" / "preprocess_summary.tsv") as f:
        f.readline()  # Headers
        stats = f.readline().split("\t")
        assert int(stats[1]) == 400  # Input reads
        assert int(stats[2]) == 400  # Both kept by cutadapt
        assert (
            int(stats[4]) + int(stats[6]) + int(stats[9]) == 200
        )  # Human + phiX + komplexity
        assert int(stats[4]) == int(stats[5])  # Human = human_copy
        assert int(stats[7]) + int(stats[9]) == 200  # Host + komplexity
        assert int(stats[8]) == 200  # Nonhost
