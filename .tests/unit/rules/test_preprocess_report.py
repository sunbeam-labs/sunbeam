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
        data_dir / "qc" / "log",
        output_dir / "sunbeam_output" / "qc" / "log",
    )
    shutil.copytree(
        data_dir / "qc" / "03_komplexity",
        output_dir / "sunbeam_output" / "qc" / "03_komplexity",
    )
    shutil.copytree(
        data_dir / "qc" / "cleaned", output_dir / "sunbeam_output" / "qc" / "cleaned"
    )
    shutil.copytree(
        data_dir / "qc" / "decontam",
        output_dir / "sunbeam_output" / "qc" / "decontam",
    )
    shutil.copytree(data_dir / "logs", output_dir / "sunbeam_output" / "logs")

    yield output_dir

    shutil.rmtree(output_dir / "sunbeam_output")


def test_preprocess_report(setup):
    output_dir = setup
    sunbeam_output_dir = output_dir / "sunbeam_output"
    r = sunbeam_output_dir / "qc" / "reports" / "preprocess_summary.tsv"

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

    assert r.stat().st_size >= 100

    with open(r) as f:
        f.readline()  # Headers
        stats = f.readline().split("\t")
        assert int(stats[1]) == 400  # Input reads
        assert int(stats[2]) == 400  # Both kept by trimmomatic
        assert (
            int(stats[6]) + int(stats[8]) + int(stats[11]) == 200
        )  # Human + phiX + komplexity
        assert int(stats[6]) == int(stats[7])  # Human = human_copy
        assert int(stats[9]) + int(stats[11]) == 200  # Host + komplexity
        assert int(stats[10]) == 200  # Nonhost
