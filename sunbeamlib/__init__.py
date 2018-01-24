__author__ = "Erik Clarke"
__license__ = "GPL2+"

import os
import re
import sys
from setuptools_scm import get_version

from semantic_version import Version
from snakemake.utils import listfiles
from snakemake.workflow import expand
from Bio import SeqIO

__version__ = str(Version.coerce(get_version()))

def build_sample_list(data_fp, filename_fmt, samplelist_fp, excluded):
    if os.path.isfile(str(samplelist_fp)):
        sys.stderr.write("Building sample list using {}...\n".format(samplelist_fp))
        Samples = _build_samples_from_file(data_fp, filename_fmt, samplelist_fp, excluded)
    else:
        sys.stderr.write("Building sample list from {}/{}...\n".format(data_fp, filename_fmt))
        Samples = _build_samples_from_dir(data_fp, filename_fmt, excluded)
    return Samples

def _build_samples_from_dir(data_fp, filename_fmt, excluded):
    """
    Build a list of samples from a data filepath and filename format.

    :param data_fp: a Path to data files
    :param filename_fmt: a string giving wildcards for {sample} and (optionally)
        {rp}, for read pair (e.g. R1 or R2).
    :param exclude: a list of sample names to exclude
    :returns: A dictionary of samples, with sample names as keys:
       Samples = {
         'sample1': {
           'R1': 'path/to/sample1_R1.fastq.gz',
           'R2': 'path/to/sample1_R2.fastq.gz',
           'paired': True
         }, ...
       }
    """
    files = list(listfiles(str(data_fp/filename_fmt)))
    Samples = {t[1]['sample']: {} for t in files if t[1]['sample'] not in excluded}
    for f in files:
        fpath = f[0]
        wcards = f[1]
        if wcards['sample'] in excluded:
            continue
        rp = wcards['rp'] if 'rp' in wcards.keys() else False
        if rp:
            Samples[wcards['sample']][rp] = fpath
            Samples[wcards['sample']]['paired'] = True
        else:
            Samples[wcards['sample']]['file'] = fpath
            Samples[wcards['sample']]['paired'] = False
    return Samples

def _build_samples_from_file(data_fp, filename_fmt, samplelist_fp, excluded):
    """
    Build a list of samples from a sample list file.
    
    :param data_fp: a Path to data files
    :param filename_fmt: a string giving wildcards for {sample} and (optionally)
        {rp}, for read pair (e.g. R1 or R2).
    :param samplelist_fp: a Path to a whitespace-delimited samplelist file,
       where the first entry is the sample name and the rest is ignored.
    :param exclude: a list of sample names to exclude
    :returns: A dictionary of samples with the same structure as
       `_build_samples_from_dir`. Paths are reconstructed using `data_fp`
       and `filename_fmt`. Samples are assumed to be paired reads.
    """
    Samples = {}
    with open(str(samplelist_fp)) as f:
        for line in f:
            sample = line.strip().split('\t')[0]
            if sample in excluded:
                continue
            Samples[sample] = {}
            for rp in ['R1', 'R2']:
                Samples[sample][rp] = _check_sample_path(
                    sample,
                    expand(str(data_fp/filename_fmt), sample=sample, rp=rp)[0])
            Samples[sample]['paired'] = True
    return Samples

def _check_sample_path(sample, fp):
    path = Path(fp)
    if not path.exists():
        sys.stderr.write(
            "Warning: original file for sample '{}' not found at {}\n".format(
                sample, fp))
    return path.resolve()

def index_files(genome, index_fp):
    """
    Return the bowtie index files for a file.
    """
    fwd, rev = (
        expand(
            "{index_fp}/{genome}.{index}.bt2",
            index_fp=index_fp,
            genome=genome,
            index=range(1,5)),
        expand(
            "{index_fp}/{genome}.rev.{index}.bt2",
            index_fp=index_fp,
            genome=genome,
            index=range(1,3))
    )
    return fwd + rev

def circular(seq, kmin, kmax, min_len):
    """Determine if a sequence is circular.

    Checks for repeated k-mer at beginning and end of a sequence for a given
    range of values for k.
    :param seq: a character sequence
    :param kmin: the smallest value of k to check
    :param kmax: the largest value of k to check
    :return: True if any overlapping k-mers found, False otherwise
    """
    if len(seq) < min_len:
        return False
    # Short-circuit checking: returns True for the first kmer that matches
    for k in range(kmin, kmax+1):
        if seq[0:k] == seq[len(seq)-k:]:
            return True
    return False

def read_seq_ids(fasta_fp):
    """
    Return the sequence identifiers for a given fasta filename.
    """
    return [record.id for record in SeqIO.parse(str(fasta_fp), "fasta")]
