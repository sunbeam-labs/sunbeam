import os
import sys
import yaml
from pathlib import Path
from pkg_resources import resource_stream
from typing import Dict, TextIO, Tuple, Union

from sunbeamlib import __version__, Version


def makepath(path: str) -> Path:
    return Path(path).expanduser()


def verify(path: str) -> Path:
    path = Path(path)
    if path.exists():
        return path.resolve()
    else:
        raise ValueError("Path %s does not exist" % path)


def validate_paths(cfg: Dict[str, str], root: Path) -> Dict[str, Union[str, Path]]:
    """Process paths in config file subsection.

    For each key ending in _fp, the value is:
    - converted to a pathlib.Path
    - ensured to be an absolute path, by appending `root` if needed
    - ensured to exist, if it is not the value from `output_fp`
    - expanded home directory ~

    :param cfg: a config file subsection
    :param root: the root directory for the project
    :returns: an updated copy of cfg
    """
    new_cfg = dict()
    for k, v in cfg.items():
        if k.endswith("_fp"):
            try:
                v = makepath(v)
            except TypeError as e:
                raise TypeError(f"Missing value for key: {k}")
            if not v.is_absolute():
                v = root / v
            if k != "output_fp":
                try:
                    v = verify(v)
                except ValueError:
                    raise ValueError("For key '%s': path '%s' does not exist" % (k, v))
        new_cfg[k] = v
    return new_cfg


def check_compatibility(cfg: Dict[str, Dict[str, str]]) -> Tuple[str, str]:
    """Returns the major version numbers from the package and config file, respectively"""

    cfg_version = Version(cfg["all"].get("version", "0.0.0"))
    pkg_version = Version(__version__)

    return (pkg_version.major, cfg_version.major)


def check_config(cfg: Dict[str, Dict[str, str]]) -> Dict[str, Dict[str, str]]:
    """Resolve root in config file, then validate paths."""

    if "root" in cfg["all"]:
        root = verify(cfg["all"]["root"])
    else:
        root = Path.cwd()
    # Iteratively check paths for each subsection
    new_cfg = dict()
    for section, values in cfg.items():
        new_cfg[section] = validate_paths(values, root)
    new_cfg["all"]["root"] = root
    return new_cfg


def output_subdir(cfg: Dict[str, Dict[str, str]], section: str) -> Path:
    return cfg["all"]["output_fp"] / cfg[section]["suffix"]


def _update_dict(
    target: Dict[str, Union[str, Dict]], new: Dict[str, Union[str, Dict]]
) -> Dict[str, Union[str, Dict]]:
    for k, v in new.items():
        if isinstance(v, dict):
            target[k] = _update_dict(target.get(k, {}), v)
        else:
            target[k] = v
    return target


def _update_dict_strict(
    target: Dict[str, Union[str, Dict]], new: Dict[str, Union[str, Dict]]
) -> Dict[str, Union[str, Dict]]:
    for k, v in new.items():
        target[k] = _update_dict_strict(target.get(k, {}), v)
    return target


def update(
    config_str: str, new: Dict[str, Union[str, Dict]], strict: bool = False
) -> Dict[str, Union[str, Dict]]:
    config = yaml.safe_load(config_str)
    if strict:
        config = _update_dict_strict(config, new)
    else:
        sbx_config = yaml.safe_load(extension_config())
        if sbx_config:
            for k, v in sbx_config.items():
                if k not in config.keys():
                    # Only update sbx_config if that sbx isn't in the config file yet
                    config = _update_dict(config, {k: v})
        config = _update_dict(config, new)
    return config


def new(
    project_fp: Union[str, Path], version: str = __version__, template: TextIO = None
) -> str:
    if template:
        config = template.read()
    else:
        config = str(
            resource_stream("sunbeamlib", "data/default_config.yml").read().decode()
        )
        # add config from extensions
        config = config + extension_config()

    return config.format(PROJECT_FP=project_fp, SB_VERSION=version)


def extension_config() -> str:
    config = ""
    sunbeam_dir = Path(os.getenv("SUNBEAM_DIR", os.getcwd()))
    for sbx in os.listdir(sunbeam_dir / "extensions"):
        if sbx[:3] != "sbx":
            continue
        try:
            sbx_files = os.listdir(sunbeam_dir / "extensions" / sbx)
        except NotADirectoryError:
            continue
        if "config.yml" in sbx_files:
            # append it to the existing config
            sbx_config_fp = sunbeam_dir / "extensions" / sbx / "config.yml"
            sbx_configfile = open(sbx_config_fp)
            sbx_config = "\n" + sbx_configfile.read()
            sbx_configfile.close()
            config = str(config + sbx_config)
    return config


def load_defaults(default_name: str) -> Dict[str, Union[str, Dict]]:
    return yaml.safe_load(
        resource_stream("sunbeamlib", "data/{}.yml".format(default_name))
        .read()
        .decode()
    )


def dump(config: Union[str, Dict], out: TextIO = sys.stdout) -> None:
    if isinstance(config, dict):
        yaml.safe_dump(config, out)
    else:
        out.write(config)
