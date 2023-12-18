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
        headers = list(map(str.strip, f.readline().split("\t")))

        stats = dict(
            zip(headers, list(map(str.strip, f.readline().split("\t"))))
        )  # LONG
        assert (
            int(stats["input"]) == 2000
        ), f"Preprocess report: Wrong number of input reads: {stats}"  # Input reads
        assert (
            int(stats["both_kept"]) == 2000
        ), f"Preprocess report: Wrong number of trimmomatic reads: {stats}"  # Both kept by trimmomatic
        assert (
            int(stats["nonhost"]) + int(stats["komplexity"]) == 2000
        ), f"Preprocess report: Wrong number of nonhost and komplexity reads: {stats}"  # Nonhost + komplexity

        stats = dict(
            zip(headers, list(map(str.strip, f.readline().split("\t"))))
        )  # SHORT
        assert (
            int(stats["input"]) == 400
        ), f"Preprocess report: Wrong number of input reads: {stats}"  # Input reads
        assert (
            int(stats["both_kept"]) == 400
        ), f"Preprocess report: Wrong number of trimmomatic reads: {stats}"  # Both kept by trimmomatic
        assert (
            int(stats["human"])
            + int(stats["phix174"])
            + int(stats["nonhost"])
            + int(stats["komplexity"])
            == 400
        ), f"Preprocess report: Wrong number of nonhost, host and komplexity reads: {stats}"  # Human + phiX + nonhost + komplexity
        assert int(stats["human"]) == int(
            stats["human_copy"]
        ), f"Preprocess report: Human doesn't equal human_copy: {stats}"  # Human = human_copy


@pytest.fixture
def init_dirty(output_dir):
    output_dir = output_dir / "sunbeam_run_dirty"

    sp.check_output(
        [
            "sunbeam",
            "init",
            "--data_fp",
            f"{test_dir / 'data' / 'dirty_reads'}",
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
            shutil.copytree(output_dir, "output_sunbeam_run_dirty/")
        except FileExistsError as e:
            pass


def test_sunbeam_run_all_dirty(init_dirty):
    output_dir = init_dirty
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

    with open(sunbeam_output_dir / "qc" / "reports" / "preprocess_summary.tsv") as f:
        f.readline()  # Headers
        stats = f.readline().split("\t")
        assert int(stats[1]) == 100  # Input reads
        assert int(stats[2]) == 22  # Both kept by trimmomatic
        assert int(stats[3]) == 0  # Fwd only
        assert int(stats[4]) == 67  # Rev only
        assert int(stats[5]) == 11  # Dropped
        assert int(stats[2]) + int(stats[3]) + int(stats[4]) + int(stats[5]) == int(
            stats[1]
        )
        assert int(stats[6]) == int(stats[7])  # Human = human_copy


@pytest.fixture
def init_no_host(output_dir):
    output_dir = output_dir / "sunbeam_run_no_host"

    sp.check_output(
        [
            "sunbeam",
            "init",
            "--data_fp",
            f"{test_dir / 'data' / 'reads'}",
            output_dir,
        ]
    )

    yield output_dir

    if os.environ.get("CI", False):
        try:
            shutil.copytree(output_dir, "output_sunbeam_run_no_host/")
        except FileExistsError as e:
            pass


def test_sunbeam_run_all_no_host(init):
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

    assert len(os.listdir(sunbeam_output_dir / "qc" / "cleaned")) == 4


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

        stats = f.readline().split("\t")  # LONG
        assert int(stats[1]) == 2000  # Input reads
        assert int(stats[2]) == 2000  # Both kept by trimmomatic
        assert int(stats[8]) + int(stats[9]) == 2000  # Nonhost + komplexity

        stats = f.readline().split("\t")  # SHORT
        assert int(stats[1]) == 400  # Input reads
        assert int(stats[2]) == 400  # Both kept by trimmomatic
        assert int(stats[4]) == int(stats[5])  # Human = human_copy
        assert int(stats[7]) == 200  # Host
        assert int(stats[8]) == 200  # Nonhost

@pytest.fixture
def init_spaces(output_dir):
    output_dir = output_dir / "sunbeam_run_space_in_header"

    sp.check_output(
        [
            "sunbeam",
            "init",
            "--data_fp",
            f"{test_dir / 'data' / 'spaces'}",
            output_dir,
        ]
    )

    yield output_dir

    if os.environ.get("CI", False):
        try:
            shutil.copytree(output_dir, "output_sunbeam_run_space_in_header/")
        except FileExistsError as e:
            pass


def test_sunbeam_run_all_space_in_header(init):
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

    assert len(os.listdir(sunbeam_output_dir / "qc" / "cleaned")) == 4
