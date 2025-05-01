import pytest
from sunbeam import CONFIGS_DIR
from sunbeam.project import SampleList, SunbeamConfig, SunbeamProfile
from sunbeam.scripts.init import main as Init
from sunbeam.scripts.run import main as Run


def test_sunbeam_run(tmp_path, DATA_DIR):
    project_dir = tmp_path / "test"
    
    Init(
        [
            str(project_dir),
            "--data_fp",
            str(DATA_DIR / "reads"),
        ]
    )

    Run(
        [
            "--profile",
            str(project_dir),
        ]
    )

    sunbeam_output = project_dir / "sunbeam_output"
    assert sunbeam_output.exists()