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
    output_decontam_dir = output_dir / "sunbeam_output" / "qc" / "decontam"
    (output_decontam_dir / "human").mkdir()
    (output_decontam_dir / "phix174").mkdir()
    shutil.copyfile(
        data_dir / "qc" / "decontam" / "human" / "unmapped_TEST.sam",
        output_decontam_dir / "human" / "unmapped_TEST.sam",
    )
    shutil.copyfile(
        data_dir / "qc" / "decontam" / "phix174" / "unmapped_TEST.sam",
        output_decontam_dir / "phix174" / "unmapped_TEST.sam",
    )

    yield output_dir

    shutil.rmtree(output_dir / "sunbeam_output")


def test_sam_convert_paired(setup):
    output_dir = setup
    sunbeam_output_dir = output_dir / "sunbeam_output"
    human1 = sunbeam_output_dir / "qc" / "decontam" / "human" / "unmapped_TEST_1.fastq"
    human2 = sunbeam_output_dir / "qc" / "decontam" / "human" / "unmapped_TEST_2.fastq"
    phix1 = sunbeam_output_dir / "qc" / "decontam" / "phix174" / "unmapped_TEST_1.fastq"
    phix2 = sunbeam_output_dir / "qc" / "decontam" / "phix174" / "unmapped_TEST_2.fastq"

    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
            "--allowed-rules=sam_convert_paired",
            "--rerun-triggers=input",
            f"{human1}",
            f"{human2}",
            f"{phix1}",
            f"{phix2}",
        ]
    )

    assert human1.stat().st_size >= 40000
    assert human2.stat().st_size >= 40000
    assert phix1.stat().st_size >= 40000
    assert phix2.stat().st_size >= 40000
