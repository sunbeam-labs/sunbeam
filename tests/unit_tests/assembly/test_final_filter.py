import shutil
import tempfile
import unittest
import os

from workflow.scripts.assembly.final_filter_f import (
    parse_fasta,
    filter_seqs,
    write_fasta,
)


class FinalFilterTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

        self.input_fp = os.path.join(self.temp_dir, "in.fa")
        self.input_contents = "> 1\nAAAAAAAAAAA\n> 2\nCCC\n> 3\nGGGGGGGG\n> 4\nTT"
        with open(self.input_fp, "w") as f:
            f.write(self.input_contents)

        self.output_fp = os.path.join(self.temp_dir, "out.fa")
        self.output_contents = "> 1\nAAAAAAAAAAA\n> 3\nGGGGGGGG\n> 2\nCCC\n"

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_all(self):
        with open(self.input_fp) as f:
            with open(self.output_fp, "w") as g:
                write_fasta(g, list(filter_seqs(parse_fasta(f), 3)))

        with open(self.output_fp) as f:
            self.assertEqual("".join(f.readlines()), self.output_contents)
