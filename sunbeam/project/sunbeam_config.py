import yaml
from pathlib import Path
from sunbeam import __version__, EXTENSIONS_DIR
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

    @classmethod
    def from_file(cls, config_fp: Path) -> "SunbeamConfig":
        """
        Create a SunbeamConfig object from a file
        """
        with open(config_fp) as f:
            config = yaml.safe_load(f)
        return cls(config)

    @classmethod
    def from_template(
        cls, template_fp: Path, root_fp: Path, extensions_dir: Path = None
    ) -> "SunbeamConfig":
        """
        Create a SunbeamConfig object from a template file
        """
        if not extensions_dir:
            extensions_dir = EXTENSIONS_DIR()

        with open(template_fp) as f:
            config = dict(yaml.safe_load(f))

        config["all"]["root"] = str(root_fp)
        config["all"]["version"] = __version__

        # Iterate over extension dirs to find config files
        for ext_name, ext_dir in cls.get_extensions(extensions_dir).items():
            config_fp = ext_dir / "config.yml"
            if config_fp.exists():
                with open(config_fp) as f:
                    ext_config = yaml.safe_load(f)
                    if ext_config:
                        config = {**config, **ext_config}

        return cls(config)

    def to_file(self, config_fp: Path):
        """
        Write the config to a file
        """
        with open(config_fp, "w") as f:
            yaml.safe_dump(self.config, f)

    def fill_missing(self, extensions_dir: Path = None):
        """
        Fill in missing extension config values with defaults
        """
        if not extensions_dir:
            extensions_dir = EXTENSIONS_DIR()

        for ext_name, ext_dir in self.get_extensions(extensions_dir).items():
            config_fp = ext_dir / "config.yml"
            if config_fp.exists():
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

        extensions = {}
        for ext_dir in extensions_dir.iterdir():
            if ext_dir.is_dir() and ext_dir.name.startswith("sbx"):
                extensions[ext_dir.name] = ext_dir

        return extensions

    @staticmethod
    def get_extension_rules(extension_fp: Path) -> list[Path]:
        """
        Find all .smk and .rules files in the extension directory using glob
        """
        return list(extension_fp.glob("**/*.smk")) + list(
            extension_fp.glob("**/*.rules")
        )

    def resolved_paths(self) -> dict[str, Path | str]:
        """
        Resolve all paths in the config file (any field ending in "_fp")
        Relative paths are resolved relative to the 'root' key
        """
        root_fp = Path(self.config["all"]["root"]).resolve()

        resolved = {}
        for key, value in self.config.items():
            if isinstance(value, dict):
                resolved[key] = {}
                for sub_key, sub_value in value.items():
                    if sub_key.endswith("_fp"):
                        if not Path(sub_value).is_absolute():
                            resolved[key][sub_key] = root_fp / sub_value
                        else:
                            resolved[key][sub_key] = Path(sub_value).resolve()
                    else:
                        resolved[key][sub_key] = sub_value
            else:
                resolved[key] = value

        resolved["all"]["root"] = root_fp
        return resolved


def output_subdir(cfg: dict[str, dict[str, str]], section: str) -> Path:
    """
    Get the output subdirectory for a given section.
    Here mostly for backwards compatibility.
    """
    return cfg["all"]["output_fp"] / section
