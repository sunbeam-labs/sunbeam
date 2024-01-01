import os
import subprocess as sp
import sys
from pathlib import Path

test_dir = Path(__file__).parent.parent.resolve()
sys.path.append(test_dir)
from config_fixture import output_dir, config


def test_sunbeam_extend():
    sp.check_output(
        [
            "sunbeam",
            "extend",
            "https://github.com/sunbeam-labs/sbx_template.git",
        ]
    )

    template_fp = Path(os.environ.get("SUNBEAM_DIR")) / "extensions" / "sbx_template"
    assert template_fp.exists()
