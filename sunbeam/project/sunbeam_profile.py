import yaml
from pathlib import Path


# Significantly limited in scope compared to SunbeamConfig
class SunbeamProfile:
    def __init__(self, config: dict = {}):
        self.config = config

    @classmethod
    def from_template(cls, template_fp: Path) -> "SunbeamProfile":
        with open(template_fp) as f:
            config = dict(yaml.safe_load(f))
        return cls(config)

    def to_file(self, config_fp: Path):
        with open(config_fp, "w") as f:
            yaml.safe_dump(self.config, f)
