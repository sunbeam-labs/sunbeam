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
def extend_template():
    sp.check_output(
        [
            "sunbeam",
            "extend",
            "https://github.com/sunbeam-labs/sbx_template.git",
        ]
    )

    yield Path(os.environ.get("SUNBEAM_DIR")) / "extensions" / "sbx_template"

    shutil.rmtree(Path(os.environ.get("SUNBEAM_DIR")) / "extensions" / "sbx_template")


def test_sunbeam_extend(extend_template):
    template_fp = extend_template

    assert template_fp.exists()
