import pytest
import subprocess as sp
import sys
from pathlib import Path
from ruamel.yaml import YAML

test_dir = Path(__file__).parent.parent.resolve()
sys.path.append(test_dir)
from config_fixture import output_dir, config


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

@pytest.fixture
def config_modify(init):
    output_dir = init
    hosts_fp = test_dir / 'data' / 'hosts'
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

    yaml = YAML(typ="safe")
    with open(output_dir / "sunbeam_config.yml") as f:
        config_dict = yaml.load(f)
    
    assert config_dict["qc"]["host_fp"] == str(hosts_fp)