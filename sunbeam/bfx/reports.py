import re
import os
import sys
from collections import OrderedDict
from io import StringIO
from typing import TextIO
import pathlib


def make_fastqc_report(
        input_report_fps,
        output_report_fp,
        log_fp,
        ):
    qual_data_long = []
    for fp in input_report_fps:
        fp = pathlib.Path(fp)
        sample_id = (fp.resolve().parent.name).split("_fastqc")[0] # Hate this
        with open(fp) as f:
            for base_idx, mean_qual in parse_fastqc_quality(f):
                qual_data_long.append((sample_id, base_idx, mean_qual))
    with open(output_report_fp, "w") as f:
        for row in simple_pivot(qual_data_long, "Samples", ""):
            f.write("\t".join(row))
            f.write("\n")


def simple_pivot(triples, id_colname, empty_val):
    # We need to iterate through triples twice
    triples = list(triples)

    colnames = {}
    colnames[id_colname] = None
    for _, colname, _ in triples:
        colnames[colname] = None
    yield list(colnames)

    nrow = 0
    current_rowname = None
    current_vals = {}
    for rowname, colname, val in triples:
        if rowname != current_rowname:
            # Emit previous row
            if nrow > 0:
                yield list(current_vals.get(k, empty_val) for k in colnames)
            # Set up current row
            current_rowname = rowname
            current_vals = {}
            current_vals[id_colname] = rowname

        current_vals[colname] = val
        nrow += 1

    # Emit final row if there were any
    if nrow > 0:
        yield list(current_vals.get(k, empty_val) for k in colnames)


def parse_trim_summary_paired(f: TextIO) -> OrderedDict[str, str]:
    for line in f.readlines():
        if line.startswith("Input Read"):
            vals = re.findall("\\D+: (\\d+)", line)
            keys = ("input", "both_kept", "fwd_only", "rev_only", "dropped")
            return OrderedDict(zip(keys, vals))


def parse_trim_summary_single(f: TextIO) -> OrderedDict[str, str]:
    for line in f:
        if line.startswith("Input Read"):
            vals = re.findall("\\D+: (\\d+)", line)
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
):
    """Return a dataframe for summary information for trimmomatic and decontam rule"""
    import pandas
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


def parse_fastqc_quality(f):
    in_per_base_quals = False
    header_seen = False
    for line in f:
        if in_per_base_quals:
            line = line.rstrip()
            if not header_seen:
                assert line.startswith("#Base")
                header_seen = True
            elif line.startswith(">>END_MODULE"):
                break
            else:
                vals = line.split("\t")
                yield (vals[0], vals[1])
        elif line.startswith(">>Per base sequence quality"):
            in_per_base_quals = True
