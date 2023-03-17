import os
import pytest
import shutil
import subprocess as sp
import sys
from pathlib import Path

test_dir = Path(__file__).parent.parent.parent.resolve()
sys.path.append(str(test_dir))
from config_fixture import output_dir, config


@pytest.fixture
def init(output_dir):
    output_dir = output_dir / "rules"
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

    (output_dir / "sunbeam_output" / "qc" / "00_samples").mkdir(
        parents=True, exist_ok=True
    )
    os.symlink(
        test_dir / "data" / "reads" / "LONG_R1.fastq.gz",
        output_dir / "sunbeam_output" / "qc" / "00_samples" / "LONG_1.fastq.gz",
    )
    os.symlink(
        test_dir / "data" / "reads" / "LONG_R2.fastq.gz",
        output_dir / "sunbeam_output" / "qc" / "00_samples" / "LONG_2.fastq.gz",
    )
    os.symlink(
        test_dir / "data" / "reads" / "SHORT_R1.fastq.gz",
        output_dir / "sunbeam_output" / "qc" / "00_samples" / "SHORT_1.fastq.gz",
    )
    os.symlink(
        test_dir / "data" / "reads" / "SHORT_R2.fastq.gz",
        output_dir / "sunbeam_output" / "qc" / "00_samples" / "SHORT_2.fastq.gz",
    )

    yield output_dir

    if os.environ.get("CI", False):
        try:
            shutil.copytree(output_dir, "output_rules/")
        except FileExistsError as e:
            pass
