import os
import pytest
import shutil
import subprocess as sp
import sys
from pathlib import Path
from . import init

test_dir = Path(__file__).parent.parent.parent.resolve()
sys.path.append(str(test_dir))
from config_fixture import output_dir, config

data_dir = Path(__file__).parent / "data"


@pytest.fixture
def setup(init):
    output_dir = init

    shutil.copytree(
        data_dir / "qc" / "01_cutadapt",
        output_dir / "sunbeam_output" / "qc" / "01_cutadapt",
    )

    yield output_dir

    shutil.rmtree(output_dir / "sunbeam_output")


def test_trimmomatic_paired(setup):
    output_dir = setup
    sunbeam_output_dir = output_dir / "sunbeam_output"
    r1 = sunbeam_output_dir / "qc" / "02_trimmomatic" / "TEST_1.fastq.gz"
    r2 = sunbeam_output_dir / "qc" / "02_trimmomatic" / "TEST_2.fastq.gz"
    ur1 = (
        sunbeam_output_dir
        / "qc"
        / "02_trimmomatic"
        / "unpaired"
        / "TEST_1_unpaired.fastq.gz"
    )
    ur2 = (
        sunbeam_output_dir
        / "qc"
        / "02_trimmomatic"
        / "unpaired"
        / "TEST_2_unpaired.fastq.gz"
    )

    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
            "--rerun-triggers=input",
            f"{r1}",
            f"{r2}",
            f"{ur1}",
            f"{ur2}",
        ]
    )

    assert r1.stat().st_size >= 10000
    assert r2.stat().st_size >= 10000
