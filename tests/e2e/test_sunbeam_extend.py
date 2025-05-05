from sunbeam import EXTENSIONS_DIR
from sunbeam.scripts.extend import main as Extend


def test_sunbeam_extend():
    Extend(["sbx_template"])

    ext_path = EXTENSIONS_DIR() / "sbx_template"
    assert ext_path.exists()


def test_sunbeam_extend_with_github_url():
    Extend(["https://github.com/sunbeam-labs/sbx_template.git"])

    ext_path = EXTENSIONS_DIR() / "sbx_template"
    assert ext_path.exists()
