__author__ = "Erik Clarke"
__license__ = "GPL2+"

import os
import re
import sys
import csv

from pathlib import Path

from semantic_version import Version
from sunbeamlib.parse import parse_fasta

__version__ = str(Version.coerce(os.environ.get("SUNBEAM_VER", "0.0.0")))


def load_sample_list(samplelist_fp, paired_end=True, root_proj=""):
    """
    Build a list of samples from a sample list file.
    :param samplelist_fp: a Path to a whitespace-delimited samplelist file,
       where the first entry is the sample name and the rest is ignored.
    :returns: A dictionary of samples with sample name and associated file(s)
    """
    Samples = {}
    with open(str(samplelist_fp)) as f:
        reader = csv.DictReader(f, fieldnames=["sample", "1", "2"])
        for row in reader:
            sample = row["sample"]
            try:
                r1 = _verify_path(row["1"])
            except ValueError:
                raise ValueError("Associated file for {} not found.".format(sample))
            r2 = ""
            if paired_end:
                try:
                    r2 = _verify_path(row["2"])
                except ValueError:
                    raise ValueError(
                        "Paired-end files specified, but mate pair for '{}' "
                        "is missing or does not exist.".format(sample)
                    )
            Samples[sample] = {"1": r1, "2": r2}
    return Samples


def guess_format_string(fnames, paired_end=True, split_pattern="([_.])"):
    """
    Try to guess the format string given a list of filenames.
    :param fnames: a list of filename strings
    :param paired_end: if True, will try to find a read-pair designator
    :param split_pattern: regex to split filenames on
    :returns: a format string with a {sample} and (maybe) an {rp} element
    """

    non_fastq_gz = [fn for fn in fnames if not fn.endswith(".fastq.gz")]
    if len(non_fastq_gz) > 0:
        raise SampleFormatError(f"Found non fastq.gz files: {str(non_fastq_gz)}")
    if len(fnames) == 0:
        raise SampleFormatError("No files in directory!")

    if not paired_end:
        return "{sample}.fastq.gz"

    splits = [list(reversed(re.split(split_pattern, fname))) for fname in fnames]
    if len(set([len(p) for p in splits])) > 1:
        sys.stderr.write(
            f"Warning: samples have inconsistent numbers of {split_pattern} characters\n"
        )

    pattern = re.compile(".+\_[1-2]\.fastq\.gz")
    if [fn for fn in fnames if pattern.match(fn)] == fnames:
        return "{sample}_{rp}.fastq.gz"
    pattern = re.compile(".+\_R[1-2]\.fastq\.gz")
    if [fn for fn in fnames if pattern.match(fn)] == fnames:
        return "{sample}_R{rp}.fastq.gz"
    pattern = re.compile("[1-2]\_.+\.fastq\.gz")
    if [fn for fn in fnames if pattern.match(fn)] == fnames:
        return "{rp}_{sample}.fastq.gz"
    pattern = re.compile("R[1-2]\_.+\.fastq\.gz")
    if [fn for fn in fnames if pattern.match(fn)] == fnames:
        return "R{rp}_{sample}.fastq.gz"

    raise SampleFormatError(
        "Couldn't identify sample format, please specify with --format option"
    )


class MissingMatePairError(Exception):
    pass


class SampleFormatError(Exception):
    pass


def _verify_path(fp):
    if not fp:
        raise ValueError("Missing filename")
    path = Path(fp)
    if not path.is_file():
        raise ValueError("File not found")
    return str(path.resolve())


def circular(seq, kmin, kmax, min_len):
    """Determine if a sequence is circular.

    Checks for repeated k-mer at beginning and end of a sequence for a given
    range of values for k."""
    if len(seq) < min_len:
        return False
    # Short-circuit checking: returns True for the first kmer that matches
    return any([k for k in range(kmin, kmax + 1) if seq[0:k] == seq[len(seq) - k :]])


def read_seq_ids(fasta_fp):
    """
    Return the sequence identifiers for a given fasta filename.
    """
    with open(str(fasta_fp)) as f:
        return list(parse_fasta(f))
