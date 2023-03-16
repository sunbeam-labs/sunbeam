import gzip
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


def test_adapter_removal_paired(setup):
    output_dir = setup
    sunbeam_output_dir = output_dir / "sunbeam_output"
    lr1 = sunbeam_output_dir / "qc" / "01_cutadapt" / "SHORT_1.fastq.gz"
    lr2 = sunbeam_output_dir / "qc" / "01_cutadapt" / "SHORT_2.fastq.gz"
    sr1 = sunbeam_output_dir / "qc" / "01_cutadapt" / "LONG_1.fastq.gz"
    sr2 = sunbeam_output_dir / "qc" / "01_cutadapt" / "LONG_2.fastq.gz"

    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
            "--allowed-rules=adapter_removal_paired",
            "--rerun-triggers=input",
            f"{lr1}",
            f"{lr2}",
            f"{sr1}",
            f"{sr2}",
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


def test_adapater_removal_paired_no_adapters(setup):
    output_dir = setup
    sunbeam_output_dir = output_dir / "sunbeam_output"
    lr1 = sunbeam_output_dir / "qc" / "01_cutadapt" / "SHORT_1.fastq.gz"
    lr2 = sunbeam_output_dir / "qc" / "01_cutadapt" / "SHORT_2.fastq.gz"
    sr1 = sunbeam_output_dir / "qc" / "01_cutadapt" / "LONG_1.fastq.gz"
    sr2 = sunbeam_output_dir / "qc" / "01_cutadapt" / "LONG_2.fastq.gz"

    config_str = f"qc: {{fwd_adapters: }}"
    sp.check_output(
        [
            "sunbeam",
            "config",
            "modify",
            "-i",
            "-s",
            config_str,
            f"{output_dir / 'sunbeam_config.yml'}",
        ]
    )
    config_str = f"qc: {{rev_adapters: }}"
    sp.check_output(
        [
            "sunbeam",
            "config",
            "modify",
            "-i",
            "-s",
            config_str,
            f"{output_dir / 'sunbeam_config.yml'}",
        ]
    )

    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
            "--allowed-rules=adapter_removal_paired",
            "--rerun-triggers=input",
            f"{lr1}",
            f"{lr2}",
            f"{sr1}",
            f"{sr2}",
        ]
    )

    assert os.path.islink(lr1)
    assert os.path.islink(lr2)
    assert os.path.islink(sr1)
    assert os.path.islink(sr2)
