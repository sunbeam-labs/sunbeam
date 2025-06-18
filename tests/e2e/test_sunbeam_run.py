import pytest
import subprocess as sp
import sys
import types
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

    ret = sp.run(
        [
            "sunbeam",
            "run",
            "--profile",
            str(project_dir),
            "--exclude",
            "all",
            "clean_qc",
            "-n",
            "--quiet",
        ],
        check=True,
        capture_output=True,
    )

    assert ret.returncode == 0
    assert "clean_qc" in ret.stderr.decode("utf-8")
    assert "filter_reads" not in ret.stderr.decode("utf-8")


def test_sunbeam_run_ai_option(tmp_path, monkeypatch):
    project_dir = tmp_path / "empty"
    project_dir.mkdir()

    called = {"flag": False}

    def dummy_create(**kwargs):
        called["flag"] = True
        return types.SimpleNamespace(
            choices=[
                types.SimpleNamespace(message=types.SimpleNamespace(content="analysis"))
            ]
        )

    fake_openai = types.SimpleNamespace()

    fake_openai.OpenAI = lambda *args, **kwargs: types.SimpleNamespace(
        chat=types.SimpleNamespace(
            completions=types.SimpleNamespace(create=dummy_create)
        )
    )

    monkeypatch.setitem(sys.modules, "openai", fake_openai)
    monkeypatch.setenv("OPENAI_API_KEY", "token")

    with pytest.raises(SystemExit):
        Run(
            [
                "--profile",
                str(project_dir),
                "--ai",
                "-n",
            ]
        )

    assert called["flag"]
