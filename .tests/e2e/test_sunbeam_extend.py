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
    output_dir = output_dir / "sunbeam_extend"

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

    yield output_dir

    if os.environ.get("CI", False):
        try:
            shutil.copytree(output_dir, "output_sunbeam_extend/")
        except FileExistsError as e:
            pass


def test_sunbeam_extend(init):
    output_dir = init
    sunbeam_output_dir = output_dir / "sunbeam_output"

    sp.check_output(
        [
            "sunbeam",
            "extend",
            "https://github.com/sunbeam-labs/sbx_template.git",
        ]
    )

    sp.check_output(
        [
            "sunbeam",
            "config",
            "update",
            "-i",
            f"{output_dir / 'sunbeam_config.yml'}",
        ]
    )

    template_fp = Path(os.environ.get("SUNBEAM_DIR")) / "extensions" / "sbx_template"
    print(template_fp)
    assert template_fp.exists()

    dag = sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            f"{output_dir}",
            "--notemp",
            "-n",
            "all_template",
        ]
    )

    assert "all_template" in dag.decode()
    assert "example_rule" in dag.decode()
    assert "example_with_script" in dag.decode()
