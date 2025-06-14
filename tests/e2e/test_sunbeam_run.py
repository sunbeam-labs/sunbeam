import pytest
from sunbeam import CONFIGS_DIR, EXTENSIONS_DIR
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

    with pytest.raises(SystemExit) as excinfo:
        Run(
            [
                "--profile",
                str(project_dir),
                "--exclude",
                "all",
            ]
        )
    assert excinfo.value.code == 0

    sunbeam_output = project_dir / "sunbeam_output"
    assert sunbeam_output.exists()


def test_sunbeam_run_with_single_end_reads(tmp_path, DATA_DIR):
    project_dir = tmp_path / "test"

    Init(
        [
            str(project_dir),
            "--data_fp",
            str(DATA_DIR / "single_end_reads"),
            "--single_end",
        ]
    )

    with pytest.raises(SystemExit) as excinfo:
        Run(
            [
                "--profile",
                str(project_dir),
                "--exclude",
                "all",
            ]
        )
    assert excinfo.value.code == 0

    sunbeam_output = project_dir / "sunbeam_output"
    assert sunbeam_output.exists()


def test_sunbeam_run_with_dirty_reads(tmp_path, DATA_DIR):
    project_dir = tmp_path / "test"

    Init(
        [
            str(project_dir),
            "--data_fp",
            str(DATA_DIR / "dirty_reads"),
        ]
    )

    with pytest.raises(SystemExit) as excinfo:
        Run(
            [
                "--profile",
                str(project_dir),
                "--exclude",
                "all",
            ]
        )
    assert excinfo.value.code == 0

    sunbeam_output = project_dir / "sunbeam_output"
    assert sunbeam_output.exists()


def test_sunbeam_run_with_target_after_exclude(tmp_path, DATA_DIR, capsys):
    project_dir = tmp_path / "test"

    Init(
        [
            str(project_dir),
            "--data_fp",
            str(DATA_DIR / "reads"),
        ]
    )

    with pytest.raises(SystemExit) as excinfo:
        Run(
            [
                "--profile",
                str(project_dir),
                "--exclude",
                "all",
                "clean_qc",
                "-n",
                "--quiet",
            ]
        )
    assert excinfo.value.code == 0

    captured = capsys.readouterr()

    assert "clean_qc" in captured.err
    assert "filter_reads" not in captured.err
