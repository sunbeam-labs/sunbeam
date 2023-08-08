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

    yield output_dir

    shutil.rmtree(output_dir / "sunbeam_output")


def test_trimmomatic_paired(setup):
    output_dir = setup
    sunbeam_output_dir = output_dir / "sunbeam_output"
    lr1 = sunbeam_output_dir / "qc" / "02_trimmomatic" / "LONG_1.fastq.gz"
    lr2 = sunbeam_output_dir / "qc" / "02_trimmomatic" / "LONG_2.fastq.gz"
    sr1 = sunbeam_output_dir / "qc" / "02_trimmomatic" / "SHORT_1.fastq.gz"
    sr2 = sunbeam_output_dir / "qc" / "02_trimmomatic" / "SHORT_2.fastq.gz"
    lur1 = (
        sunbeam_output_dir
        / "qc"
        / "02_trimmomatic"
        / "unpaired"
        / "LONG_1_unpaired.fastq.gz"
    )
    lur2 = (
        sunbeam_output_dir
        / "qc"
        / "02_trimmomatic"
        / "unpaired"
        / "LONG_2_unpaired.fastq.gz"
    )
    sur1 = (
        sunbeam_output_dir
        / "qc"
        / "02_trimmomatic"
        / "unpaired"
        / "SHORT_1_unpaired.fastq.gz"
    )
    sur2 = (
        sunbeam_output_dir
        / "qc"
        / "02_trimmomatic"
        / "unpaired"
        / "SHORT_2_unpaired.fastq.gz"
    )

    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
            "--allowed-rules=trimmomatic_paired",
            "--rerun-triggers=input",
            f"{lr1}",
            f"{lr2}",
            f"{lur1}",
            f"{lur2}",
            f"{sr1}",
            f"{sr2}",
            f"{sur1}",
            f"{sur2}",
        ]
    )

    assert lr1.stat().st_size >= 50000
    assert lr2.stat().st_size >= 50000
    assert sr1.stat().st_size >= 10000
    assert sr2.stat().st_size >= 10000

    with gzip.open(lr1) as f1, gzip.open(lr2) as f2:
        assert len(f1.readlines()) == len(f2.readlines())
    with gzip.open(sr1) as f1, gzip.open(sr2) as f2:
        assert len(f1.readlines()) == len(f2.readlines())
