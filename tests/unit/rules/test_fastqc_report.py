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
        data_dir / "qc" / "reports", output_dir / "sunbeam_output" / "qc" / "reports"
    )
    os.remove(output_dir / "sunbeam_output" / "qc" / "reports" / "fastqc_quality.tsv")

    yield output_dir

    shutil.rmtree(output_dir / "sunbeam_output")


def test_fastqc_report(setup):
    output_dir = setup
    sunbeam_output_dir = output_dir / "sunbeam_output"
    r = sunbeam_output_dir / "qc" / "reports" / "fastqc_quality.tsv"

    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
            "--allowed-rules=fastqc_report",
            "--rerun-triggers=input",
            f"{r}",
        ]
    )

    with open(r) as f:
        assert True
