import os
import pytest
import shutil
import sys
import tempfile
import yaml
from pathlib import Path


@pytest.fixture(autouse=True)
def DATA_DIR():
    return Path(__file__).parent / "data"


@pytest.fixture
def config():
    tests_dir = Path(__file__).parent.resolve()

    with open(tests_dir / "test_config.yml") as f:
        config_dict = yaml.safe_load(f)

    yield config_dict


@pytest.fixture
def output_dir(config):
    config = config
    output_dir = Path()

    if not config["output_dir"]:
        output_dir = Path(tempfile.mkdtemp())
    else:
        output_dir = Path(config["output_dir"])
        output_dir.mkdir(parents=True, exist_ok=True)
        if not os.listdir(output_dir) == []:
            if config["overwrite"]:
                shutil.rmtree(output_dir)
                output_dir.mkdir()
            else:
                sys.exit(
                    "overwrite is set to false but output_dir points to a non-empty directory"
                )

    if not config["temp_env"]:
        pass
    else:
        # TODO: Create temp_env
        pass

    yield output_dir

    if not config["output_dir"]:
        shutil.rmtree(output_dir)
