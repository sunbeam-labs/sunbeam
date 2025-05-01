import os
import pytest
import shutil
import subprocess as sp
import sys
import yaml
from pathlib import Path

test_dir = Path(__file__).parent.parent.resolve()
sys.path.append(test_dir)
from tests.conftest import output_dir, config


@pytest.fixture
def init(output_dir):
    output_dir = output_dir / "sunbeam_config"

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
            shutil.copytree(output_dir, "output_sunbeam_config/")
        except FileExistsError as e:
            pass


@pytest.fixture
def config_modify(init):
    output_dir = init
    hosts_fp = test_dir / "data" / "hosts"
    config_str = f"qc: {{host_fp: {hosts_fp}}}"

    sp.check_output(
        [
            "sunbeam",
            "config",
            "modify",
            "-i",
            "-s",
            f"{config_str}",
            f"{output_dir / 'sunbeam_config.yml'}",
        ]
    )

    yield output_dir, hosts_fp


def test_config_modify(config_modify):
    output_dir, hosts_fp = config_modify

    with open(output_dir / "sunbeam_config.yml") as f:
        config_dict = yaml.safe_load(f)

    assert config_dict["qc"]["host_fp"] == str(hosts_fp)


@pytest.fixture
def config_update(init):
    output_dir = init
