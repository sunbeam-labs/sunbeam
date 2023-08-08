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
        data_dir / "qc" / "cleaned", output_dir / "sunbeam_output" / "qc" / "cleaned"
    )

    yield output_dir

    shutil.rmtree(output_dir / "sunbeam_output")


def test_align_to_host(setup):
    output_dir = setup
    sunbeam_output_dir = output_dir / "sunbeam_output"
    lhuman = (
        sunbeam_output_dir / "qc" / "decontam" / "intermediates" / "human" / "LONG.sam"
    )
    lphix = (
        sunbeam_output_dir
        / "qc"
        / "decontam"
        / "intermediates"
        / "phix174"
        / "LONG.sam"
    )
    shuman = (
        sunbeam_output_dir / "qc" / "decontam" / "intermediates" / "human" / "SHORT.sam"
    )
    sphix = (
        sunbeam_output_dir
        / "qc"
        / "decontam"
        / "intermediates"
        / "phix174"
        / "SHORT.sam"
    )

    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
            "--allowed-rules=align_to_host",
            "--rerun-triggers=input",
            f"{lhuman}",
            f"{lphix}",
            f"{shuman}",
            f"{sphix}",
        ]
    )

    assert lhuman.stat().st_size >= 500000
    assert lphix.stat().st_size >= 500000
    assert shuman.stat().st_size >= 100000
    assert sphix.stat().st_size >= 100000
