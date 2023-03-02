import os
import pytest
import shutil
import sys
from pathlib import Path

test_dir = Path(__file__).parent.parent.parent.resolve()
sys.path.append(str(test_dir))

from config_fixture import output_dir, config
from sunbeamlib import (
    load_sample_list,
    guess_format_string,
    _verify_path,
    circular,
    read_seq_ids,
)

data_dir = Path(__file__).parent / "data"


@pytest.fixture
def init(output_dir):
    output_dir = output_dir / "sunbeamlib"
    output_dir.mkdir(parents=True, exist_ok=True)

    yield output_dir

    if os.environ.get("CI", False):
        try:
            shutil.copytree(output_dir, "output_sunbeamlib/")
        except FileExistsError as e:
            pass


def test_load_sample_list(init):
    output_dir = init


def test_guess_format_string(init):
    output_dir = init

    test_strs = ["TEST_R1.fastq.gz", "TEST_R2.fastq.gz"]
    ret = guess_format_string(test_strs)
    assert ret == "{sample}_R{rp}.fastq.gz"
