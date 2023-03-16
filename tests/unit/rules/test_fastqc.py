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
    lr1 = sunbeam_output_dir / "qc" / "reports" / "LONG_1_fastqc" / "fastqc_data.txt"
    lr2 = sunbeam_output_dir / "qc" / "reports" / "LONG_2_fastqc" / "fastqc_data.txt"
    sr1 = sunbeam_output_dir / "qc" / "reports" / "SHORT_1_fastqc" / "fastqc_data.txt"
    sr2 = sunbeam_output_dir / "qc" / "reports" / "SHORT_2_fastqc" / "fastqc_data.txt"

    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
            "--allowed-rules=fastqc",
            "--rerun-triggers=input",
            f"{lr1}",
            f"{lr2}",
            f"{sr1}",
            f"{sr2}",
        ]
    )

    assert lr1.stat().st_size >= 5000
    assert lr2.stat().st_size >= 5000
    assert sr1.stat().st_size >= 30000
    assert sr2.stat().st_size >= 30000
