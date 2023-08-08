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
        data_dir / "qc" / "decontam" / "intermediates",
        output_dir / "sunbeam_output" / "qc" / "decontam" / "intermediates",
    )
    shutil.copytree(
        data_dir / "qc" / "cleaned", output_dir / "sunbeam_output" / "qc" / "cleaned"
    )

    yield output_dir

    shutil.rmtree(output_dir / "sunbeam_output")


def test_filter_reads(setup):
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

    clr1 = sunbeam_output_dir / "qc" / "cleaned" / "LONG_1.fastq.gz"
    clr2 = sunbeam_output_dir / "qc" / "cleaned" / "LONG_2.fastq.gz"
    csr1 = sunbeam_output_dir / "qc" / "cleaned" / "SHORT_1.fastq.gz"
    csr2 = sunbeam_output_dir / "qc" / "cleaned" / "SHORT_2.fastq.gz"

    out = sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
            "--allowed-rules=filter_reads",
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

    assert clr1.stat().st_size == lr1.stat().st_size
    assert clr2.stat().st_size == lr2.stat().st_size

    with open(ll1) as f:
        assert f.readline() == "human\thuman_copy\tphix174\thost\tnonhost\n"
        assert f.readline() == "0\t0\t0\t0\t1995\n"

    with gzip.open(lr1) as f1, gzip.open(lr2) as f2:
        assert len(f1.readlines()) == len(f2.readlines())

    assert sr1.stat().st_size >= 5000
    assert sr2.stat().st_size >= 5000

    assert csr1.stat().st_size > sr1.stat().st_size
    assert csr2.stat().st_size > sr2.stat().st_size

    with open(sl1) as f:
        assert f.readline() == "human\thuman_copy\tphix174\thost\tnonhost\n"
        assert f.readline() == "94\t94\t100\t194\t198\n"

    with gzip.open(sr1) as f1, gzip.open(sr2) as f2:
        assert len(f1.readlines()) == len(f2.readlines())
