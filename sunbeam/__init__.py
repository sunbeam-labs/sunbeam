__author__ = "Erik Clarke"
__license__ = "GPL2+"

from snakemake.utils import listfiles
from snakemake.workflow import expand


def build_sample_list(data_fp, filename_fmt, excluded):
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


def index_files(genome):
    """
    Return the bowtie index files for a file.
    """
    fwd, rev = (
        expand(
            "{index_fp}/{genome}.{index}.bt2",
            index_fp=Cfg['bt2_index_fp'],
            genome=genome,
            index=range(1,5)),
        expand(
            "{index_fp}/{genome}.rev.{index}.bt2",
            index_fp=Cfg['bt2_index_fp'],
            genome=genome,
            index=range(1,3))
    )
    return fwd + rev
