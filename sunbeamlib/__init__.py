__author__ = "Erik Clarke"
__license__ = "GPL2+"

from snakemake.utils import listfiles
from snakemake.workflow import expand
import os


def build_sample_list(data_fp, filename_fmt, excluded):
    if os.path.isdir(str(data_fp)):
        Samples = build_sample_from_dir(data_fp, filename_fmt, excluded)
    else:
        Samples = build_sample_from_file(data_fp)
    return Samples

def build_sample_from_dir(data_fp, filename_fmt, excluded):
    """
    Build a list of samples from a data filepath and filename format.

    :param data_fp: a Path to data files
    :param filename_fmt: a string giving wildcards for {sample} and (optionally)
        {rp}, for read pair (e.g. R1 or R2).
    :param exclude: a list of sample names to exclude
    :returns: A dictionary of samples, with sample names as keys
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

def build_sample_from_file(data_fp):
    """"
    Build a list of samples from a barcode file
    :param bc_file: a Path to barcode file
    :returns: A dictionary of samples, with sample names as keys
    """
    with open(str(data_fp)) as f:
        lines = f.read().splitlines()
    ids = []
    for line in lines:
         ids.append(line.split("\t")[0])
    # todo: not sure about adding the path of actual reads
    Samples = dict((id,"paired") for id in ids)
    return Samples

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

