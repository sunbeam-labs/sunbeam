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
        data_dir / "qc" / "decontam" / "human" / "unmapped_LONG.sam",
        output_decontam_dir / "human" / "unmapped_LONG.sam",
    )
    shutil.copyfile(
        data_dir / "qc" / "decontam" / "human" / "unmapped_SHORT.sam",
        output_decontam_dir / "human" / "unmapped_SHORT.sam",
    )
    shutil.copyfile(
        data_dir / "qc" / "decontam" / "phix174" / "unmapped_LONG.sam",
        output_decontam_dir / "phix174" / "unmapped_LONG.sam",
    )
    shutil.copyfile(
        data_dir / "qc" / "decontam" / "phix174" / "unmapped_SHORT.sam",
        output_decontam_dir / "phix174" / "unmapped_SHORT.sam",
    )

    yield output_dir

    shutil.rmtree(output_dir / "sunbeam_output")


def test_sam_convert_paired(setup):
    output_dir = setup
    sunbeam_output_dir = output_dir / "sunbeam_output"
    lhuman1 = sunbeam_output_dir / "qc" / "decontam" / "human" / "unmapped_LONG_1.fastq"
    lhuman2 = sunbeam_output_dir / "qc" / "decontam" / "human" / "unmapped_LONG_2.fastq"
    lphix1 = sunbeam_output_dir / "qc" / "decontam" / "phix174" / "unmapped_LONG_1.fastq"
    lphix2 = sunbeam_output_dir / "qc" / "decontam" / "phix174" / "unmapped_LONG_2.fastq"
    shuman1 = sunbeam_output_dir / "qc" / "decontam" / "human" / "unmapped_SHORT_1.fastq"
    shuman2 = sunbeam_output_dir / "qc" / "decontam" / "human" / "unmapped_SHORT_2.fastq"
    sphix1 = sunbeam_output_dir / "qc" / "decontam" / "phix174" / "unmapped_SHORT_1.fastq"
    sphix2 = sunbeam_output_dir / "qc" / "decontam" / "phix174" / "unmapped_SHORT_2.fastq"

    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
            "--allowed-rules=sam_convert_paired",
            "--rerun-triggers=input",
            f"{lhuman1}",
            f"{lhuman2}",
            f"{lphix1}",
            f"{lphix2}",
            f"{shuman1}",
            f"{shuman2}",
            f"{sphix1}",
            f"{sphix2}",
        ]
    )

    assert lhuman1.stat().st_size >= 300000
    assert lhuman2.stat().st_size >= 300000
    assert lphix1.stat().st_size >= 300000
    assert lphix2.stat().st_size >= 300000
    assert shuman1.stat().st_size >= 40000
    assert shuman2.stat().st_size >= 40000
    assert sphix1.stat().st_size >= 40000
    assert sphix2.stat().st_size >= 40000
