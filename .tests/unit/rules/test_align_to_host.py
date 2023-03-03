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
    shutil.copytree(
        data_dir / "qc" / "02_trimmomatic",
        output_dir / "sunbeam_output" / "qc" / "02_trimmomatic",
    )
    shutil.copytree(
        data_dir / "qc" / "log" / "komplexity",
        output_dir / "sunbeam_output" / "qc" / "log" / "komplexity",
    )
    shutil.copytree(
        data_dir / "qc" / "03_komplexity",
        output_dir / "sunbeam_output" / "qc" / "03_komplexity",
    )
    shutil.copytree(
        data_dir / "qc" / "cleaned", output_dir / "sunbeam_output" / "qc" / "cleaned"
    )

    yield output_dir

    shutil.rmtree(output_dir / "sunbeam_output")


def test_align_to_host(setup):
    output_dir = setup
    sunbeam_output_dir = output_dir / "sunbeam_output"
    human = (
        sunbeam_output_dir / "qc" / "decontam" / "intermediates" / "human" / "TEST.sam"
    )
    phix = (
        sunbeam_output_dir
        / "qc"
        / "decontam"
        / "intermediates"
        / "phix174"
        / "TEST.sam"
    )

    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
            "--rerun-triggers=input",
            f"{human}",
            f"{phix}",
        ]
    )

    assert human.stat().st_size >= 100000
    assert phix.stat().st_size >= 100000
