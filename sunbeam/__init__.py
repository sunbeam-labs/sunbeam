import os
from pathlib import Path


__version__ = "5.0.0.c"
__author__ = "Erik Clarke"
__license__ = "GPL2+"


EXTENSIONS_DIR = lambda: Path(
    os.environ.get("SUNBEAM_EXTENSIONS", Path(__file__).parent.resolve() / "extensions")
)
WORKFLOW_DIR = Path(__file__).parent.resolve() / "workflow"
CONFIGS_DIR = Path(__file__).parent.resolve() / "configs"


def get_docker_str(repo: str, user: str = "sunbeamlabs") -> str:
    docker_tag = os.environ.get("SUNBEAM_DOCKER_TAG", f"v{__version__}")

    return f"docker://{user}/{repo}:{docker_tag}"


def get_ext_path(ext_name: str) -> Path:
    ext_path = EXTENSIONS_DIR() / ext_name

    if ext_path.exists():
        return ext_path
    raise ValueError(f"Extension {ext_name} not found in {EXTENSIONS_DIR()}.")
