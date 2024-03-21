import os
import pytest
import shutil
import sys
from pathlib import Path

test_dir = Path(__file__).parent.parent.parent.resolve()
sys.path.append(str(test_dir))

from config_fixture import output_dir, config
from sunbeamlib.config import *

data_dir = Path(__file__).parent / "data"


@pytest.fixture
def init(output_dir):
    output_dir = output_dir / "sunbeamlib"
    output_dir.mkdir(parents=True, exist_ok=True)

    yield output_dir

    if os.environ.get("CI", False):
        try:
            shutil.copytree(output_dir, "output_sunbeamlib/")
        except FileExistsError as e:
            pass


def test_makepath_on_none():
    assert makepath(None) == Path("")
    assert makepath("test") == Path("test")