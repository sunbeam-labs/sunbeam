import shutil
import tempfile
import unittest
import os

from scripts.assembly.get_coverage_f import parse_depth, get_cov_stats, write_csv


class GetCoverageTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

        self.input_fp = os.path.join(self.temp_dir, "in.depth")
        self.input_contents = "k141_1\t1\t3\nk141_1\t2\t3\nk141_1\t3\t3\nk141_0\t1\t1\nk141_0\t2\t3\nk141_0\t3\t3"
        with open(self.input_fp, "w") as f:
            f.write(self.input_contents)

        self.output_fp = os.path.join(self.temp_dir, "out.csv")
        self.output_contents = """sample,contig,sum,min,max,mean,median,stddev,coverage,length
ecoli,k141_1,9,3,3,3.0,3.0,0.0,3,3
ecoli,k141_0,7,1,3,2.3333333333333335,3.0,0.9428090415820634,3,3\n"""

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_all(self):
        with open(self.input_fp) as f:
            with open(self.output_fp, "w") as g:
                write_csv(g, get_cov_stats(parse_depth(f), "ecoli"))

        with open(self.output_fp) as f:
            self.assertEqual("".join(f.readlines()), self.output_contents)
