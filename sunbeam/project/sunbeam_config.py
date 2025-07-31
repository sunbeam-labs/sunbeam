import os
import yaml
from pathlib import Path
from sunbeam import __version__, EXTENSIONS_DIR, logger
from typing import Dict, Union


class SunbeamConfig:
    """
    Sunbeam configuration file

    Defining samples:
      Run 'sunbeam list_samples <data_dir>' to create a list of samples and
      associated fastq.gz files. Samples must be in gzipped fastq format.

    Paths:
      Paths are resolved through the following rules:
        1. If the path is absolute, the path is parsed as-is
        2. If the path is not absolute, the path at 'root' is appended to it
        3. If the path is not 'output_fp', the path is checked to ensure it exists

    Suffixes:
      Each subsection contains a 'suffix' key that defines the folder under
       'output_fp' where the results of that section are put.
    """

    def __init__(self, config: Dict[str, Union[str, Dict]] = {}):
        self.config = config
        logger.debug(f"Initialized SunbeamConfig with keys: {list(self.config.keys())}")

    @classmethod
    def from_file(cls, config_fp: Path) -> "SunbeamConfig":
        """
        Create a SunbeamConfig object from a file
        """
        logger.debug(f"Loading config from file: {config_fp}")
        with open(config_fp) as f:
            config = yaml.safe_load(f)

        if "all" in config and "root" in config["all"]:
            root = Path(config["all"]["root"]).resolve()
            config["all"]["root"] = str(root)

        logger.debug(f"Loaded config keys: {list(config.keys())}")
        return cls(config)

    @classmethod
    def from_template(
        cls, template_fp: Path, root_fp: Path, extensions_dir: Path = None
    ) -> "SunbeamConfig":
        """
        Create a SunbeamConfig object from a template file
        """
        logger.debug(f"Creating config from template {template_fp} with root {root_fp}")
        if not extensions_dir:
            extensions_dir = EXTENSIONS_DIR()

        with open(template_fp) as f:
            config = dict(yaml.safe_load(f))

        config["all"]["root"] = str(Path(root_fp).resolve())
        config["all"]["version"] = __version__

        # Iterate over extension dirs to find config files
        for ext_name, ext_dir in cls.get_extensions(extensions_dir).items():
            config_fp = ext_dir / "config.yml"
            if config_fp.exists():
                logger.debug(f"Loading extension config from {config_fp}")
                with open(config_fp) as f:
                    ext_config = yaml.safe_load(f)
                    if ext_config:
                        config = {**config, **ext_config}

        logger.debug("Config from template loaded")
        return cls(config)

    def to_file(self, config_fp: Path):
        """
        Write the config to a file
        """
        logger.debug(f"Writing config to {config_fp}")
        with open(config_fp, "w") as f:
            yaml.safe_dump(self.config, f)

    def fill_missing(self, extensions_dir: Path = None):
        """
        Fill in missing extension config values with defaults
        """
        logger.debug(
            f"Filling missing config values using extensions dir: {extensions_dir}"
        )
        if not extensions_dir:
            extensions_dir = EXTENSIONS_DIR()

        for ext_name, ext_dir in self.get_extensions(extensions_dir).items():
            config_fp = ext_dir / "config.yml"
            if config_fp.exists():
                logger.debug(f"Checking defaults in {config_fp}")
                with open(config_fp) as f:
                    ext_config = yaml.safe_load(f)
                    if ext_config:
                        # Accounting for the possibility of multiple configs in a single file
                        for key, value in ext_config.items():
                            if key not in self.config:
                                self.config[key] = value
                            else:
                                for sub_key, sub_value in value.items():
                                    if sub_key not in self.config[key]:
                                        self.config[key][sub_key] = sub_value

    def modify(self, change_str: str):
        """
        Modify the config file with the specified changes
        change_str should be a string in the format "root_key: {sub_key: value}"
        """
        logger.debug(f"Applying config modification: {change_str}")
        changes = yaml.safe_load(change_str)
        for k, v in changes.items():
            if k not in self.config:
                self.config[k] = v
            else:
                if isinstance(v, dict):
                    for sub_k, sub_v in v.items():
                        self.config[k][sub_k] = sub_v
                else:
                    self.config[k] = v

    @staticmethod
    def get_extensions(extensions_dir: Path = None) -> dict[str, Path]:
        """
        Get a list of all extensions in the extensions directory
        """
        if not extensions_dir:
            extensions_dir = EXTENSIONS_DIR()
        logger.debug(f"Scanning for extensions in {extensions_dir}")

        extensions = {}
        for ext_dir in extensions_dir.iterdir():
            if ext_dir.is_dir() and ext_dir.name.startswith("sbx"):
                extensions[ext_dir.name] = ext_dir
                logger.debug(f"Found extension: {ext_dir.name}")

        return extensions

    @staticmethod
    def get_extension_rules(extension_fp: Path) -> list[Path]:
        """
        Find all .smk and .rules files in the extension directory using glob
        """
        logger.debug(f"Searching for rule files in {extension_fp}")
        return list(extension_fp.glob("**/*.smk")) + list(
            extension_fp.glob("**/*.rules")
        )

    def resolved_paths(self) -> dict[str, Path | str]:
        """
        Resolve all paths in the config file (any field ending in "_fp")
        Relative paths are resolved relative to the 'root' key
        """
        logger.debug("Resolving configuration paths")
        root_fp = Path(self.config["all"]["root"])
        logger.debug(f"Root path resolved to {root_fp}")

        resolved = {}
        for key, value in self.config.items():
            if isinstance(value, dict):
                resolved[key] = {}
                for sub_key, sub_value in value.items():
                    if sub_key.endswith("_fp"):
                        if not sub_value:
                            logger.warning(
                                f"{key}.{sub_key} is empty, setting to os.devnull"
                            )
                            resolved[key][sub_key] = Path(os.devnull)
                        elif not Path(sub_value).is_absolute():
                            resolved[key][sub_key] = root_fp / sub_value
                            logger.debug(
                                f"Resolved {key}.{sub_key} to {resolved[key][sub_key]}"
                            )
                        else:
                            resolved[key][sub_key] = Path(sub_value).resolve()
                            logger.debug(
                                f"Resolved {key}.{sub_key} to {resolved[key][sub_key]}"
                            )
                    else:
                        resolved[key][sub_key] = sub_value
            else:
                resolved[key] = value

        resolved["all"]["root"] = root_fp
        logger.debug("All paths resolved")
        return resolved


def output_subdir(cfg: dict[str, dict[str, str]], section: str) -> Path:
    """
    Get the output subdirectory for a given section.
    Here mostly for backwards compatibility.
    """
    subdir = cfg["all"]["output_fp"] / section
    logger.debug(f"Output subdir for {section}: {subdir}")
    return subdir
