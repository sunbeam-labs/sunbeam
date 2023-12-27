import pandas
import re
import os
import sys
from collections import OrderedDict
from io import StringIO
from typing import TextIO


def parse_trim_summary_paired(f: TextIO) -> OrderedDict[str, str]:
    for line in f.readlines():
        if line.startswith("Input Read"):
            vals = re.findall("\D+\: (\d+)", line)
            keys = ("input", "both_kept", "fwd_only", "rev_only", "dropped")
            return OrderedDict(zip(keys, vals))


def parse_trim_summary_single(f: TextIO) -> OrderedDict[str, str]:
    for line in f:
        if line.startswith("Input Read"):
            vals = re.findall("\D+\: (\d+)", line)
            keys = ("input", "kept", "dropped")
            return OrderedDict(zip(keys, vals))


def parse_decontam_log(f: TextIO) -> OrderedDict[str, str]:
    keys = f.readline().rstrip().split("\t")
    vals = f.readline().rstrip().split("\t")
    return OrderedDict(zip(keys, vals))


def parse_komplexity_log(f: TextIO) -> OrderedDict[str, str]:
    return OrderedDict([("komplexity", str(len(f.readlines())))])


def summarize_qual_decontam(
    tfile: str, dfile: str, kfile: str, paired_end: bool
) -> pandas.DataFrame:
    """Return a dataframe for summary information for trimmomatic and decontam rule"""
    tname = os.path.basename(tfile).split(".out")[0]
    with open(tfile) as tf:
        with open(dfile) as jf:
            with open(kfile) as kf:
                if paired_end:
                    trim_data = parse_trim_summary_paired(tf)
                else:
                    trim_data = parse_trim_summary_single(tf)

                decontam_data = parse_decontam_log(jf)

                komplexity_data = parse_komplexity_log(kf)
    sys.stderr.write("trim data: {}\n".format(trim_data))
    sys.stderr.write("decontam data: {}\n".format(decontam_data))
    sys.stderr.write("komplexity data: {}\n".format(komplexity_data))
    return pandas.DataFrame(
        OrderedDict(trim_data, **(decontam_data), **(komplexity_data)),
        index=[tname.replace("trimmomatic_", "").replace(".log", "")],
    )


def parse_fastqc_quality(filename: str) -> pandas.DataFrame:
    with open(filename) as f:
        report = f.read()
    try:
        tableString = re.search(
            "\>\>Per base sequence quality.*?\n(.*?)\n\>\>END_MODULE", report, re.DOTALL
        ).group(1)

        f_s = StringIO(tableString)
        df = pandas.read_csv(
            f_s, sep="\t", usecols=["#Base", "Mean"], index_col="#Base"
        )
        sample_name = os.path.basename(filename.split("_fastqc")[0])
        df.columns = [sample_name]
        f_s.close()
        return df
    except AttributeError as e:
        sys.stderr.write(f"{filename} has no per-base sequence quality reports.")
        return None
