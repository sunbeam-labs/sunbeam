import os
import pytest
import shutil
import subprocess as sp
import sys
from pathlib import Path

test_dir = Path(__file__).parent.parent.resolve()
sys.path.append(test_dir)
from config_fixture import output_dir, config


@pytest.fixture
def init(output_dir):
    output_dir = output_dir / "sunbeam_init"

    sp.check_output(
        [
            "sunbeam",
            "init",
            "--data_fp",
            f"{test_dir / 'data' / 'reads'}",
            output_dir,
        ]
    )

    yield output_dir

    if os.environ.get("CI", False):
        try:
            shutil.copytree(output_dir, "output_sunbeam_init/")
        except FileExistsError as e:
            pass


def test_init(init):
    output_dir = init

    config_fp = output_dir / "sunbeam_config.yml"
    profile_fp = output_dir / "config.yaml"
    samples_fp = output_dir / "samples.csv"

    assert config_fp.exists()
    assert profile_fp.exists()
    assert samples_fp.exists()
