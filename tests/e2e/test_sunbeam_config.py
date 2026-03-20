from sunbeam import CONFIGS_DIR, EXTENSIONS_DIR
from sunbeam.project import SunbeamConfig
from sunbeam.scripts import Config


def test_sunbeam_config_update(tmp_path):
    project_fp = tmp_path / "test"
    project_fp.mkdir(parents=True, exist_ok=True)
    config_fp = project_fp / "sunbeam_config.yml"

    sc = SunbeamConfig.from_template(CONFIGS_DIR / "default_config.yml", project_fp)
    sc.to_file(config_fp)

    ext_fp = EXTENSIONS_DIR() / "sbx_test"
    ext_fp.mkdir(parents=True, exist_ok=True)
    with open(ext_fp / "config.yml", "w") as f:
        f.write(
            """
            sbx_test:
                test: "test"
            """
        )

    Config(
        [
            str(config_fp),
            "--update",
        ]
    )

    sc = SunbeamConfig.from_file(config_fp)
    assert sc.config["sbx_test"]["test"] == "test"


def test_sunbeam_config_modify(tmp_path):
    project_fp = tmp_path / "test"
    project_fp.mkdir(parents=True, exist_ok=True)
    config_fp = project_fp / "sunbeam_config.yml"

    sc = SunbeamConfig.from_template(CONFIGS_DIR / "default_config.yml", project_fp)
    sc.to_file(config_fp)

    Config(
        [
            str(config_fp),
            "--modify",
            "sbx_test: {test: 'modified'}",
        ]
    )

    sc = SunbeamConfig.from_file(config_fp)
    assert sc.config["sbx_test"]["test"] == "modified"
