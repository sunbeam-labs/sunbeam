import gzip
import os
import pytest
import shutil
import sys
from pathlib import Path

test_dir = Path(__file__).parent.parent.parent.resolve()
sys.path.append(str(test_dir))

from config_fixture import output_dir, config
from sunbeamlib import parse

data_dir = Path(__file__).parent.parent.parent / "data"


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
    
def test_parse_fasta():
    with open(data_dir / "hosts" / "phix174.fasta") as f:
        for lines in parse.parse_fasta(f):
            assert len(lines) == 2
            assert lines[0] == "gi|9626372|ref|NC_001422.1| Enterobacteria phage phiX174 sensu lato, complete genome"
            assert len(lines[1]) > 100
            assert set(list(lines[1])) == set(["A", "C", "G", "T"])

def test_parse_fastq():
    with gzip.open(data_dir / "reads" / "TEST_R1.fastq.gz") as f:
        for lines in parse.parse_fastq(f):
            assert len(lines) == 4
            assert set(list(lines[1])).issubset(set(["A", "C", "G", "T"]))
            assert lines[2] == "+"



