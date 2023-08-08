import os
import pytest
import shutil
import sys
import tempfile
from pathlib import Path
from ruamel.yaml import YAML


@pytest.fixture
def config():
    tests_dir = Path(__file__).parent.resolve()

    yaml = YAML(typ="safe")
    with open(tests_dir / "test_config.yml") as f:
        config_dict = yaml.load(f)

    if not shutil.which("mamba"):
        sys.exit(
            "Sunbeam needs mamba to be installed to run tests `conda install -c conda-forge mamba`. If you want to run sunbeam (not the tests) without mamba you can use `sunbeam run --profile /path/to/project/ --conda-frontend=conda`."
        )

    yield config_dict


@pytest.fixture
def output_dir(config):
    yaml = config
    output_dir = Path()

    if not yaml["output_dir"]:
        output_dir = Path(tempfile.mkdtemp())
    else:
        output_dir = Path(yaml["output_dir"])
        output_dir.mkdir(parents=True, exist_ok=True)
        if not os.listdir(output_dir) == []:
            if yaml["overwrite"]:
                shutil.rmtree(output_dir)
                output_dir.mkdir()
            else:
                sys.exit(
                    "overwrite is set to false but output_dir points to a non-empty directory"
                )

    if not yaml["temp_env"]:
        pass
    else:
        # TODO: Create temp_env
        pass

    sunbeam_dir = Path(os.environ.get("SUNBEAM_DIR"))
    extensions_fp = sunbeam_dir / "extensions"
    extensions_moved_fp = sunbeam_dir / "extensions_moved"
    try:
        shutil.move(extensions_fp, extensions_moved_fp)
        extensions_fp.mkdir()
        with open(extensions_fp / ".test", "w") as f:
            f.write("")
    except FileNotFoundError as e:
        pass

    yield output_dir

    if (extensions_fp / ".test").exists():
        shutil.rmtree(extensions_fp)
    try:
        shutil.move(extensions_moved_fp, extensions_fp)
    except FileNotFoundError as e:
        pass
    if not yaml["output_dir"]:
        shutil.rmtree(output_dir)
