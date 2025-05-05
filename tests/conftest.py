import pytest
import shutil
from pathlib import Path


@pytest.fixture(autouse=True)
def DATA_DIR() -> Path:
    return Path(__file__).parent / "data"


@pytest.fixture(autouse=True)
def set_extension_dir(monkeypatch, tmp_path):
    ext_dir = tmp_path / "extensions"
    ext_dir.mkdir(parents=True, exist_ok=True)
    monkeypatch.setenv("SUNBEAM_EXTENSIONS", str(ext_dir))
    print("Patching SUNBEAM_EXTENSIONS to", str(ext_dir))

    yield ext_dir

    shutil.rmtree(ext_dir)
    monkeypatch.undo()
