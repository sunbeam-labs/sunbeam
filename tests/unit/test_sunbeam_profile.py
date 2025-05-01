import yaml
from sunbeam import CONFIGS_DIR
from sunbeam.project.sunbeam_profile import SunbeamProfile



def test_empty_profile():
    sp = SunbeamProfile()
    assert sp.config == {}


def test_profile_from_template():
    template_fp = CONFIGS_DIR / "default_profile.yaml"
    sp = SunbeamProfile.from_template(template_fp)
    
    assert "default-resources" in sp.config


def test_profile_to_file(tmp_path):
    config_fp = tmp_path / "test_profile.yaml"
    sp = SunbeamProfile()
    sp.config = {"test-key": "test-value"}
    sp.to_file(config_fp)
    
    with open(config_fp) as f:
        loaded_config = yaml.safe_load(f)
    
    assert loaded_config == {"test-key": "test-value"}