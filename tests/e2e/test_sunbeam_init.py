from sunbeam.project.sample_list import SampleList
from sunbeam.project.sunbeam_config import SunbeamConfig
from sunbeam.scripts.init import main


def test_sunbeam_init(tmp_path):
    project_dir = tmp_path / "test"

    main(
        [
            str(project_dir),
        ]
    )

    assert (project_dir / "sunbeam_config.yml").exists()
    assert (project_dir / "config.yaml").exists()


def test_sunbeam_init_with_data(tmp_path, DATA_DIR):
    project_dir = tmp_path / "test"

    main(
        [
            str(project_dir),
            "--data_fp",
            str(DATA_DIR / "reads"),
        ]
    )

    assert (project_dir / "sunbeam_config.yml").exists()
    assert (project_dir / "config.yaml").exists()
    assert (project_dir / "samples.csv").exists()

    sl = SampleList(project_dir / "samples.csv")
    assert sl.paired_end
    assert len(sl.samples) == 2


def test_sunbeam_init_with_single_end(tmp_path, DATA_DIR):
    project_dir = tmp_path / "test"

    main(
        [
            str(project_dir),
            "--data_fp",
            str(DATA_DIR / "single_end_reads"),
            "--single_end",
        ]
    )

    assert (project_dir / "sunbeam_config.yml").exists()
    assert (project_dir / "config.yaml").exists()
    assert (project_dir / "samples.csv").exists()

    sl = SampleList(project_dir / "samples.csv", paired_end=False)
    assert len(sl.samples) == 2


def test_sunbeam_init_with_format(tmp_path, DATA_DIR):
    project_dir = tmp_path / "test"

    main(
        [
            str(project_dir),
            "--data_fp",
            str(DATA_DIR / "reads"),
            "--format",
            "{sample}_R{rp}.fastq.gz",
        ]
    )

    assert (project_dir / "sunbeam_config.yml").exists()
    assert (project_dir / "config.yaml").exists()
    assert (project_dir / "samples.csv").exists()

    sl = SampleList(project_dir / "samples.csv")
    assert len(sl.samples) == 2


def test_sunbeam_init_with_template(tmp_path):
    project_dir = tmp_path / "test"
    template_fp = tmp_path / "template_config.yml"
    with open(template_fp, "w") as f:
        f.write(
            """
            all:
                root: "{PROJECT_FP}"
                version: 0.1
            """
        )

    main(
        [
            str(project_dir),
            "--template",
            str(template_fp),
        ]
    )

    sc = SunbeamConfig.from_file(project_dir / "sunbeam_config.yml")
    assert sc.config["all"]["root"] == str(project_dir)