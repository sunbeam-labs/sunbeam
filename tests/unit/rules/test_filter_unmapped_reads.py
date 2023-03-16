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
        data_dir / "qc" / "cleaned", output_dir / "sunbeam_output" / "qc" / "cleaned"
    )
    shutil.copytree(
        data_dir / "qc" / "decontam" / "human",
        output_dir / "sunbeam_output" / "qc" / "decontam" / "human",
    )
    shutil.copytree(
        data_dir / "qc" / "decontam" / "human_copy",
        output_dir / "sunbeam_output" / "qc" / "decontam" / "human_copy",
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
    lr1 = sunbeam_output_dir / "qc" / "decontam" / "LONG_1.fastq.gz"
    lr2 = sunbeam_output_dir / "qc" / "decontam" / "LONG_2.fastq.gz"
    ll1 = sunbeam_output_dir / "qc" / "log" / "decontam" / "LONG_1.txt"
    ll2 = sunbeam_output_dir / "qc" / "log" / "decontam" / "LONG_2.txt"
    sr1 = sunbeam_output_dir / "qc" / "decontam" / "SHORT_1.fastq.gz"
    sr2 = sunbeam_output_dir / "qc" / "decontam" / "SHORT_2.fastq.gz"
    sl1 = sunbeam_output_dir / "qc" / "log" / "decontam" / "SHORT_1.txt"
    sl2 = sunbeam_output_dir / "qc" / "log" / "decontam" / "SHORT_2.txt"

    out = sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
            "--allowed-rules=filter_unmapped_reads",
            "--rerun-triggers=input",
            f"{lr1}",
            f"{lr2}",
            f"{ll1}",
            f"{ll2}",
            f"{sr1}",
            f"{sr2}",
            f"{sl1}",
            f"{sl2}",
        ]
    )

    assert lr1.stat().st_size >= 50000
    assert lr2.stat().st_size >= 50000
    assert ll1.stat().st_size >= 10
    assert ll2.stat().st_size >= 10

    with gzip.open(lr1) as f1, gzip.open(lr2) as f2:
        assert len(f1.readlines()) == len(f2.readlines())
    
    assert sr1.stat().st_size >= 5000
    assert sr2.stat().st_size >= 5000
    assert sl1.stat().st_size >= 10
    assert sl2.stat().st_size >= 10

    with gzip.open(sr1) as f1, gzip.open(sr2) as f2:
        assert len(f1.readlines()) == len(f2.readlines())
