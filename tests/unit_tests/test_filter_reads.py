import shutil
import tempfile
import unittest
import os

from collections import OrderedDict

from workflow.scripts.filter_reads_f import count_host_reads, calculate_counts, write_log


class FinalFilterTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

        self.input1_fp = os.path.join(self.temp_dir, "in1.ids")
        self.input1_contents = "read0\nread2"
        with open(self.input1_fp, "w") as f:
            f.write(self.input1_contents)

        self.input2_fp = os.path.join(self.temp_dir, "in2.ids")
        self.input2_contents = "read0\nread4"
        with open(self.input2_fp, "w") as f:
            f.write(self.input2_contents)

        self.reads_fp = os.path.join(self.temp_dir, "reads.fa")
        self.reads_contents = """@read0\nAAAAAA\n+\n++++++
@read1\nAAAAAA\n+\n++++++
@read2\nAAAAAA\n+\n++++++
@read3\nAAAAAA\n+\n++++++
@read4\nAAAAAA\n+\n++++++"""
        with open(self.reads_fp, "w") as f:
            f.write(self.reads_contents)
        os.system(f"gzip {self.reads_fp}")

        self.output_fp = os.path.join(self.temp_dir, "out.log")
        self.output_contents = "1\t3\t1\n"

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_all(self):
        with open(self.output_fp, "w") as f:
            hostdict = OrderedDict()
            net_hostlist = set()

            count_host_reads(self.input1_fp, hostdict, net_hostlist)
            count_host_reads(self.input2_fp, hostdict, net_hostlist)

            host, nonhost = calculate_counts(self.reads_fp, net_hostlist)

            write_log(f, hostdict, host, nonhost)

        with open(self.output_fp) as f:
            f.readline()
            self.assertEqual(f.readline(), self.output_contents)
