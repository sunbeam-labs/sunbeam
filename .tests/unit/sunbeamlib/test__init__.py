import os
import pytest
import shutil
import sys
from pathlib import Path

test_dir = Path(__file__).parent.parent.parent.resolve()
sys.path.append(str(test_dir))

from config_fixture import output_dir, config
from sunbeamlib import (
    __version__,
    load_sample_list,
    guess_format_string,
    _verify_path,
    circular,
    read_seq_ids,
)

data_dir = Path(__file__).parent / "data"


@pytest.fixture
def init(output_dir):
    output_dir = output_dir / "sunbeamlib"
    output_dir.mkdir(parents=True, exist_ok=True)

    yield output_dir

    if os.environ.get("CI", False):
        try:
            shutil.copytree(output_dir, "output_sunbeamlib/")
        except FileExistsError as e:
            pass

def test_version():
    assert __version__ != "0.0.0"

def test_load_sample_list(init):
    output_dir = init
    samples_fp = output_dir / "samples"
    samples_fp.mkdir()
    sample1 = samples_fp / "TEST_R1.fastq.gz"
    sample2 = samples_fp / "TEST_R2.fastq.gz"
    sample_list_fp = output_dir / "samples.csv"
    
    with open(sample_list_fp, "w") as f:
        f.write(f"TEST,{sample1.resolve()},{sample2.resolve()}")
    
    try:
        load_sample_list(sample_list_fp)
        assert False
    except ValueError as e:
        pass

    with open(sample1, "w") as f1:
        f1.write(" ")

    try:
        load_sample_list(sample_list_fp)
        assert False
    except ValueError as e:
        pass
    
    assert load_sample_list(sample_list_fp, False)["TEST"] == {"1": str(sample1.resolve()), "2": ""}

    with open(sample2, "w") as f2:
        f2.write(" ")
    
    sample_list = load_sample_list(sample_list_fp)
    assert sample_list["TEST"] == {"1": str(sample1.resolve()), "2": str(sample2.resolve())}


def test_guess_format_string():
    test_strs = ["TEST_R1.fastq.gz", "TEST_R2.fastq.gz"]
    ret = guess_format_string(test_strs)
    assert ret == "{sample}_R{rp}.fastq.gz"

def test_guess_format_string_single_end():
    test_strs = ["TEST.fastq.gz", "Test.fastq.gz", "test.fastq.gz"]
    ret = guess_format_string(test_strs, False)
    assert ret == "{sample}{sample}.fastq.gz"

def test_verify_path(init):
    output_dir = init

    try:
        _verify_path("")
        assert False
    except ValueError as e:
        pass

    try:
        _verify_path("thisdoesnotexist")
        assert False
    except ValueError as e:
        pass

    with open(output_dir / "test", "w") as f:
        f.write(" ")
    assert _verify_path(output_dir / "test") == str((output_dir / "test").resolve())
    os.remove(output_dir / "test")

def test_circular():
    assert not circular("aattccgg", 1, 3, 3)
    assert circular("aaattccggaaa", 3, 5, 2)
    assert not circular("aaattccggaaa", 4, 7, 3)