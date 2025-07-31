import os
import pytest
import yaml
from pathlib import Path
from sunbeam import __version__, CONFIGS_DIR
from sunbeam.project.sunbeam_config import SunbeamConfig


@pytest.fixture
def test_extension(tmp_path) -> Path:
    ext_name = "sbx_test_extension"
    ext_dir = tmp_path / ext_name
    ext_dir.mkdir(parents=True, exist_ok=True)
    ext_config_fp = ext_dir / "config.yml"
    ext_config = {
        "bloop": "blorp",
        "sbx_test_extension": {
            "bloop": "blorp",
            "blap": "blip",
        },
    }
    with open(ext_config_fp, "w") as f:
        yaml.dump(ext_config, f)

    return ext_dir


def test_empty_config():
    sc = SunbeamConfig()
    assert sc.config == {}


def test_config_from_template(tmp_path):
    template_fp = CONFIGS_DIR / "default_config.yml"
    root_fp = Path("root") / "root"
    ext_dir = tmp_path / "extensions"
    ext_dir.mkdir(parents=True, exist_ok=True)

    sc = SunbeamConfig.from_template(template_fp, root_fp, ext_dir)
    assert sc.config["all"]["root"] == str(root_fp.resolve())
    assert sc.config["all"]["version"] == __version__


def test_config_with_extension(test_extension):
    template_fp = CONFIGS_DIR / "default_config.yml"
    root_fp = Path("root") / "root"
    ext_dir = test_extension

    sc = SunbeamConfig.from_template(template_fp, root_fp, ext_dir.parent)

    assert len(sc.config["sbx_test_extension"]) == 2
    assert sc.config["sbx_test_extension"]["bloop"] == "blorp"


def test_config_from_file(tmp_path):
    config_str = """
    all:
        root: /path/to/root
        version: 1.0.0
    """
    config_fp = tmp_path / "config.yml"
    with open(config_fp, "w") as f:
        f.write(config_str)

    sc = SunbeamConfig.from_file(config_fp)
    assert sc.config["all"]["root"] == "/path/to/root"


def test_config_to_file(tmp_path):
    sc = SunbeamConfig(
        {
            "all": {
                "root": "/path/to/root",
                "version": __version__,
            }
        }
    )

    config_fp = tmp_path / "config.yml"
    sc.to_file(config_fp)

    with open(config_fp) as f:
        config = yaml.safe_load(f)
        assert config["all"]["root"] == "/path/to/root"


def test_fill_missing(test_extension):
    root_fp = Path("root") / "root"
    ext_dir = test_extension

    sc = SunbeamConfig(
        {
            "all": {
                "root": str(root_fp),
                "version": __version__,
            }
        }
    )

    assert "sbx_test_extension" not in sc.config

    sc.fill_missing(ext_dir.parent)

    assert len(sc.config["sbx_test_extension"]) == 2


def test_fill_missing_doesnt_overwrite(test_extension):
    root_fp = Path("root") / "root"
    ext_dir = test_extension

    sc = SunbeamConfig(
        {
            "all": {
                "root": str(root_fp),
                "version": __version__,
            },
            "sbx_test_extension": {
                "bloop": "not_blorp",
            },
        }
    )

    assert sc.config["sbx_test_extension"]["bloop"] == "not_blorp"

    sc.fill_missing(ext_dir.parent)

    assert sc.config["sbx_test_extension"]["bloop"] == "not_blorp"


def test_modify():
    config = {
        "all": {
            "root": "/path/to/root",
            "version": __version__,
        },
        "sbx_test_extension": {
            "bloop": "blorp",
        },
    }
    sc = SunbeamConfig(config)

    sc.modify("sbx_test_extension: {bloop: 'new_value'}")
    assert sc.config["sbx_test_extension"]["bloop"] == "new_value"


def test_modify_one_of_many():
    config = {
        "all": {
            "root": "/path/to/root",
            "version": __version__,
        },
        "sbx_test_extension": {
            "bloop": "blorp",
            "blap": "blip",
        },
    }
    sc = SunbeamConfig(config)

    sc.modify("sbx_test_extension: {bloop: 'new_value'}")
    assert sc.config["sbx_test_extension"]["bloop"] == "new_value"
    assert sc.config["sbx_test_extension"]["blap"] == "blip"


def test_empty_fp(tmp_path):
    config_fp = tmp_path / "config.yml"

    with open(config_fp, "w") as f:
        f.write(
            """
            all:
                root: /path/to/root
                version: 1.0.0
            sbx_test_extension:
                bloop_fp: /path/to/bloop_fp
                blap_fp: ""
                blip_fp: 
        """
        )

    sc = SunbeamConfig.from_file(config_fp)

    assert sc.config["sbx_test_extension"]["bloop_fp"] == "/path/to/bloop_fp"
    assert sc.config["sbx_test_extension"]["blap_fp"] == ""
    assert sc.config["sbx_test_extension"]["blip_fp"] is None

    res = sc.resolved_paths()

    assert res["sbx_test_extension"]["bloop_fp"] == Path("/path/to/bloop_fp")
    assert res["sbx_test_extension"]["blap_fp"] == Path(os.devnull)
    assert res["sbx_test_extension"]["blip_fp"] == Path(os.devnull)
