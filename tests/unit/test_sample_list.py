import pytest
from pathlib import Path
from sunbeam.project.sample_list import SampleList


def test_load_from_dir(DATA_DIR):
    sl = SampleList(DATA_DIR / "reads")

    assert sl.paired_end
    assert Path(sl.samples["LONG"]["r1"]).exists()

    sl = SampleList(DATA_DIR / "dirty_reads")
    sl = SampleList(DATA_DIR / "empty")
    sl = SampleList(DATA_DIR / "single_end_reads", paired_end=False)
    sl = SampleList(DATA_DIR / "spaces")


def test_to_file(tmp_path, DATA_DIR):
    sl = SampleList(DATA_DIR / "reads")
    output_fp = tmp_path / "sample_list.csv"
    sl.to_file(output_fp)

    with open(output_fp, "r") as f:
        lines = f.readlines()

    assert len(lines) == 2
    assert lines[1].strip().startswith("LONG") or lines[1].strip().startswith("SHORT")


def test_bad_sample_names(tmp_path):
    sl = SampleList()
    output_fp = tmp_path / "sample_list.csv"
    sl.samples["BAD SAMPLE NAME"] = {"r1": "test_r1.fastq.gz", "r2": "test_r2.fastq.gz"}
    sl.samples["bad_sample/name"] = {
        "r1": "test2_r1.fastq.gz",
        "r2": "test2_r2.fastq.gz",
    }
    sl.to_file(output_fp)

    with pytest.raises(ValueError):
        SampleList(output_fp)


def test_duplicate_samples(tmp_path):
    sl = SampleList()
    output_fp = tmp_path / "sample_list.csv"
    sl.samples["LONG"] = {"r1": "test_r1.fastq.gz", "r2": "test_r2.fastq.gz"}
    sl.to_file(output_fp)

    with open(output_fp, "a") as f:
        f.write("LONG,test1_r1.fastq.gz,test1_r2.fastq.gz\n")

    sl = SampleList(output_fp)
    assert len(sl.samples) == 1
