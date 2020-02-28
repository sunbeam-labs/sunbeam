import os
from unittest import TestCase

from sunbeamlib.mapping import FlagstatsParser


class TestFlagstatsParser(TestCase):

    def test_parse(self):
        file_path = os.path.join(os.path.dirname(__file__), "test_file.flagstats")
        expected_list = [
            "6013502", "0",
            "7562", "0",
            "0", "0",
            "0", "0",
            "5991648", "0",
            "6005940", "0",
            "3002970", "0",
            "3002970", "0",
            "5911116", "0",
            "5978452", "0",
            "5634", "0",
            "65898", "0",
            "36723", "0",
        ]
        tested_list = FlagstatsParser(file_path).get_list()
        self.assertListEqual(tested_list, expected_list)
