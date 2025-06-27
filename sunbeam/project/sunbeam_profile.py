import yaml
from pathlib import Path
from sunbeam import logger


# Significantly limited in scope compared to SunbeamConfig
class SunbeamProfile:
    def __init__(self, config: dict = {}):
        self.config = config
        logger.debug(
            f"Initialized SunbeamProfile with keys: {list(self.config.keys())}"
        )

    @classmethod
    def from_template(cls, template_fp: Path) -> "SunbeamProfile":
        logger.debug(f"Loading profile template from {template_fp}")
        with open(template_fp) as f:
            config = dict(yaml.safe_load(f))
        logger.debug("Profile template loaded")
        return cls(config)

    def to_file(self, config_fp: Path):
        logger.debug(f"Writing profile to {config_fp}")
        with open(config_fp, "w") as f:
            yaml.safe_dump(self.config, f)
