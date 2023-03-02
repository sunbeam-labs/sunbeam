import pytest
import subprocess as sp
import sys
from pathlib import Path

test_dir = Path(__file__).parent.parent.resolve()
sys.path.append(test_dir)
from config_fixture import output_dir, config


@pytest.fixture
def list_samples(output_dir):
    output_dir = output_dir / "sunbeam_list_samples"

    yield sp.check_output(["sunbeam", "list_samples", f"{test_dir / 'data' / 'reads'}"])


def test_list_samples(list_samples):
    sample_list = list_samples

    assert "tests/data/reads/TEST_R1.fastq.gz" in sample_list.decode("utf-8")
    assert "tests/data/reads/TEST_R2.fastq.gz" in sample_list.decode("utf-8")
