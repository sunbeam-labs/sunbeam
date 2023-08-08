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

    yield output_dir

    shutil.rmtree(output_dir / "sunbeam_output")


def test_sample_intake(init):
    output_dir = init
    sunbeam_output_dir = output_dir / "sunbeam_output"
    lr1 = sunbeam_output_dir / "qc" / "00_samples" / "LONG_1.fastq.gz"
    lr2 = sunbeam_output_dir / "qc" / "00_samples" / "LONG_2.fastq.gz"
    sr1 = sunbeam_output_dir / "qc" / "00_samples" / "SHORT_1.fastq.gz"
    sr2 = sunbeam_output_dir / "qc" / "00_samples" / "SHORT_2.fastq.gz"

    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
            f"{lr1}",
            f"{lr2}",
            f"{sr1}",
            f"{sr2}",
        ]
    )

    assert os.path.islink(lr1)
    assert os.path.islink(lr2)
    assert os.path.islink(sr1)
    assert os.path.islink(sr2)
