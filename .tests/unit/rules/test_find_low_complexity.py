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

    yield output_dir

    shutil.rmtree(output_dir / "sunbeam_output")


def test_find_low_complexity(setup):
    output_dir = setup
    sunbeam_output_dir = output_dir / "sunbeam_output"
    r = sunbeam_output_dir / "qc" / "log" / "komplexity" / "TEST.filtered_ids"

    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
            "--rerun-triggers=input",
            f"{r}",
        ]
    )

    expected_ids = [
        "KI270762.1_45667_46245_1:0:0_0:0:0_6/1",
        "KI270762.1_10179_10639_7:0:0_2:0:0_22/1",
        "KI270762.1_1730_2260_0:1:0_0:0:0_31/1",
        "KI270762.1_10398_10889_1:0:0_2:0:0_3b/1",
        "KI270762.1_13312_13765_3:0:0_0:0:0_60/1",
        "KI270762.1_1651_2168_0:0:0_0:0:0_36/2",
    ]

    with open(r) as f:
        ids = [id.strip() for id in f.readlines()]
        for id in expected_ids:
            assert id in ids
