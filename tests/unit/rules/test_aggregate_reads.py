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


def test_aggregate_reads(setup):
    output_dir = setup
    sunbeam_output_dir = output_dir / "sunbeam_output"
    shostreads = (
        sunbeam_output_dir / "qc" / "decontam" / "intermediates" / "SHORT_hostreads.ids"
    )
    lhostreads = (
        sunbeam_output_dir / "qc" / "decontam" / "intermediates" / "LONG_hostreads.ids"
    )

    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
            "--allowed-rules=aggregate_reads",
            "--rerun-triggers=input",
            f"{shostreads}",
            f"{lhostreads}",
        ]
    )

    assert shostreads.stat().st_size >= 300000
    assert lhostreads.stat().st_size >= 300000
