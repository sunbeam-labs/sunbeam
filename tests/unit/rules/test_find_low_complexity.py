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


def test_find_low_complexity(setup):
    output_dir = setup
    sunbeam_output_dir = output_dir / "sunbeam_output"
    lr = sunbeam_output_dir / "qc" / "log" / "komplexity" / "LONG.filtered_ids"
    sr = sunbeam_output_dir / "qc" / "log" / "komplexity" / "SHORT.filtered_ids"

    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
            "--allowed-rules=find_low_complexity",
            "--rerun-triggers=input",
            f"{lr}",
            f"{sr}",
        ]
    )

    expected_long_ids = [
        "NZ_CP069563.1_58419_58958_1:0:0_0:0:0_3a1/1",
        "NZ_CP069563.1_58420_58979_0:0:0_0:0:0_3c3/1",
        "NC_000913.3_59044_59553_0:1:0_0:0:0_1f0/1",
        "NC_000913.3_41948_42399_1:0:0_1:0:0_257/1",
        "NC_000913.3_58612_59129_2:0:0_2:1:0_2b/2",
    ]

    expected_short_ids = [
        "NZ_CP069563.1_58552_59002_1:0:0_1:0:0_0/2",
        "NZ_CP069563.1_30995_31446_1:0:0_0:0:0_4e/2",
        "KI270762.1_36773_37232_0:0:0_1:0:0_6/2",
        "KI270762.1_3610_4174_0:0:0_0:0:0_7/2",
        "KI270762.1_10301_10852_0:0:0_2:0:0_1f/2",
        "KI270762.1_54637_55107_1:0:0_0:0:0_22/2",
        "KI270762.1_35081_35616_3:0:0_1:0:0_2d/2",
        "KI270762.1_45718_46259_2:0:0_0:0:0_4e/2",
    ]

    with open(lr) as f:
        ids = [id.strip() for id in f.readlines()]
        for id in expected_long_ids:
            assert id in ids

    with open(sr) as f:
        ids = [id.strip() for id in f.readlines()]
        for id in expected_short_ids:
            assert id in ids
