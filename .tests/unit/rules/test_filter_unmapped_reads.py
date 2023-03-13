import gzip
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
    shutil.copytree(
        data_dir / "qc" / "decontam" / "intermediates",
        output_dir / "sunbeam_output" / "qc" / "decontam" / "intermediates",
    )
    shutil.copytree(
        data_dir / "qc" / "decontam" / "human",
        output_dir / "sunbeam_output" / "qc" / "decontam" / "human",
    )
    shutil.copytree(
        data_dir / "qc" / "decontam" / "phix174",
        output_dir / "sunbeam_output" / "qc" / "decontam" / "phix174",
    )

    yield output_dir

    shutil.rmtree(output_dir / "sunbeam_output")


def test_filter_unmapped_reads(setup):
    output_dir = setup
    sunbeam_output_dir = output_dir / "sunbeam_output"
    r1 = sunbeam_output_dir / "qc" / "decontam" / "TEST_1.fastq.gz"
    r2 = sunbeam_output_dir / "qc" / "decontam" / "TEST_2.fastq.gz"
    l1 = sunbeam_output_dir / "qc" / "log" / "decontam" / "TEST_1.txt"
    l2 = sunbeam_output_dir / "qc" / "log" / "decontam" / "TEST_2.txt"

    out = sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
            "--allowed-rules=filter_unmapped_reads",
            "--rerun-triggers=code",
            f"{r1}",
            f"{r2}",
            f"{l1}",
            f"{l2}",
            "-n"
        ]
    )

    assert r1.stat().st_size >= 5000
    assert r2.stat().st_size >= 5000
    assert l1.stat().st_size >= 10
    assert l2.stat().st_size >= 10

    with gzip.open(r1) as f1, gzip.open(r2) as f2:
        assert len(f1.readlines()) == len(f2.readlines())
