import os
from pathlib import Path


__version__ = "4.7.0"
__author__ = "Erik Clarke"
__license__ = "GPL2+"


EXTENSIONS_DIR = Path(os.environ.get("SUNBEAM_EXTENSIONS", Path(__file__).parent.resolve() / "extensions"))
WORKFLOW_DIR = Path(__file__).parent.resolve() / "workflow"
CONFIGS_DIR = Path(__file__).parent.resolve() / "configs"


def get_docker_str(repo: str, user: str = "sunbeamlabs") -> str:
    docker_tag = os.environ.get("SUNBEAM_DOCKER_TAG", f"v{__version__}")

    return f"docker://{user}/{repo}:{docker_tag}"


def get_ext_path(ext_name: str) -> Path:
    try:
        ext_path = Path(os.environ["SUNBEAM_DIR"]) / "extensions" / ext_name
    except KeyError:
        raise ValueError("SUNBEAM_DIR not set in environment.")

    if ext_path.exists():
        return ext_path
    raise ValueError(f"Extension {ext_name} not found in {os.environ['SUNBEAM_DIR']}.")
