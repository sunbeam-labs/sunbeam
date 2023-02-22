import pathlib
import pytest
from ruamel.yaml import YAML


@pytest.fixture
def config():
    tests_dir = pathlib.Path(__file__).parent.resolve()

    yaml = YAML(typ="safe")
    with open(tests_dir / "test_config.yml") as f:
        config_dict = yaml.load(f)

    yield config_dict
