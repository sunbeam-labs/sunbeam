# Home to project management code (config file handling, sample gathering, etc.)
from pathlib import Path
from sunbeam.project.sample_list import SampleList
from sunbeam.project.sunbeam_config import SunbeamConfig
from sunbeam.project.sunbeam_profile import SunbeamProfile


def output_subdir(cfg: dict[str, dict[str, str]], section: str) -> Path:
    return cfg["all"]["output_fp"] / cfg[section]["suffix"]
