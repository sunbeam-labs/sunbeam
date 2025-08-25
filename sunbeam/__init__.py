import os
from pathlib import Path
from sunbeam.logging import get_sunbeam_logger


__author__ = "Erik Clarke"
__license__ = "GPL2+"
logger = get_sunbeam_logger()


# Directory where extensions are stored. This defaults to ``~/.sunbeam/extensions``
# but can be overridden by setting the ``SUNBEAM_EXTENSIONS`` environment
# variable.  The directory is created if it does not already exist so that
# subsequent operations (such as cloning an extension) do not fail.
def EXTENSIONS_DIR() -> Path:
    default_dir = Path.home() / ".sunbeam" / "extensions"
    ext_dir = Path(os.environ.get("SUNBEAM_EXTENSIONS", default_dir))
    ext_dir.mkdir(parents=True, exist_ok=True)
    return ext_dir
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


def get_ext_version(ext_name: str) -> str:
    ext_path = get_ext_path(ext_name)
    version_file = ext_path / "VERSION"

    if version_file.exists():
        with open(version_file, "r") as f:
            return f.read().strip()
    else:
        print("Version file not found for extension:", ext_name)
        return "0.0.0"


from . import _version

__version__ = _version.get_versions()["version"]
