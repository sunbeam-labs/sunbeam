import os
import pytest
import shutil
import subprocess as sp
import sys
import tempfile
from pathlib import Path

test_dir = Path(__file__).parent.parent.resolve()
sys.path.append(test_dir)
from config_fixture import config


@pytest.fixture
def init(config):
    yaml = config
    output_dir = Path()

    if not yaml["output_dir"]:
        output_dir = Path(tempfile.mkdtemp())
    else:
        output_dir = Path(yaml["output_dir"])
        output_dir.mkdir(parents=True, exist_ok=True)
        if not os.listdir(output_dir) == 0:
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

    sp.check_output(
        [
            "sunbeam",
            "init",
            "--data_fp",
            f"{test_dir / 'data' / 'reads'}",
            output_dir,
        ]
    )

    config_str = f"qc: {{host_fp: {test_dir / 'data' / 'hosts'}}}"

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

    sunbeam_dir = Path(os.environ.get("SUNBEAM_DIR"))
    shutil.move(sunbeam_dir / "extensions", sunbeam_dir / "extensions_moved")

    yield output_dir

    shutil.move(sunbeam_dir / "extensions_moved", sunbeam_dir / "extensions")
    shutil.rmtree(output_dir)


def test_sunbeam_all(init):
    output_dir = init

    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
        ]
    )

    with open(test_dir / "targets.txt") as f:
        for line in f.readlines():
            if not line.strip():
                continue
            target = output_dir / line.strip()
            if not target.exists():
                raise SystemExit(f"Target '{target}' not found")
            elif target.stat().st_size == 0:
                raise SystemExit(f"Target '{target}' is empty")
            else:
                print(f"Found target '{target}'")
