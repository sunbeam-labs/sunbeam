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
        data_dir / "qc" / "decontam" / "intermediates",
        output_dir / "sunbeam_output" / "qc" / "decontam" / "intermediates",
    )

    yield output_dir

    shutil.rmtree(output_dir / "sunbeam_output")


def test_get_unmapped_reads(setup):
    output_dir = setup
    sunbeam_output_dir = output_dir / "sunbeam_output"
    lhuman = sunbeam_output_dir / "qc" / "decontam" / "human" / "unmapped_LONG.sam"
    lphix = sunbeam_output_dir / "qc" / "decontam" / "phix174" / "unmapped_LONG.sam"
    shuman = sunbeam_output_dir / "qc" / "decontam" / "human" / "unmapped_SHORT.sam"
    sphix = sunbeam_output_dir / "qc" / "decontam" / "phix174" / "unmapped_SHORT.sam"

    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
            "--allowed-rules=get_unmapped_reads",
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
