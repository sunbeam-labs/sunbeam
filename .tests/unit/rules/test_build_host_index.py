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


def test_build_host_index(setup):
    output_dir = setup
    sunbeam_output_dir = output_dir / "sunbeam_output"
    exts = ["amb", "ann", "bwt", "pac", "sa"]
    human = [test_dir / "data" / "hosts" / f"human.fasta.{ext}" for ext in exts]
    human_copy = [
        test_dir / "data" / "hosts" / f"human_copy.fasta.{ext}" for ext in exts
    ]
    phix = [test_dir / "data" / "hosts" / f"phix174.fasta.{ext}" for ext in exts]

    args = [
        "sunbeam",
        "run",
        "--profile",
        f"{output_dir}",
        "--notemp",
        "--allowed-rules=build_host_index",
        "--rerun-triggers=input",
    ]
    args += human
    args += human_copy
    args += phix
    sp.check_output(args)

    assert len(os.listdir(test_dir / "data" / "hosts")) == 18
