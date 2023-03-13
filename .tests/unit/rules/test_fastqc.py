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
        data_dir / "qc" / "02_trimmomatic",
        output_dir / "sunbeam_output" / "qc" / "02_trimmomatic",
    )

    yield output_dir

    shutil.rmtree(output_dir / "sunbeam_output")


def test_fastqc(setup):
    output_dir = setup
    sunbeam_output_dir = output_dir / "sunbeam_output"
    r1 = sunbeam_output_dir / "qc" / "reports" / "TEST_1_fastqc" / "fastqc_data.txt"
    r2 = sunbeam_output_dir / "qc" / "reports" / "TEST_2_fastqc" / "fastqc_data.txt"

    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
            "--allowed-rules=fastqc",
            "--rerun-triggers=input",
            f"{r1}",
            f"{r2}",
        ]
    )

    assert r1.stat().st_size >= 30000
    assert r2.stat().st_size >= 30000
