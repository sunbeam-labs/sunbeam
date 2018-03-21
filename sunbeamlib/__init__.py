__author__ = "Erik Clarke"
__license__ = "GPL2+"

import os
import re
import sys
import csv

from pathlib import Path
from pkg_resources import get_distribution

from semantic_version import Version
from snakemake.utils import listfiles
from snakemake.workflow import expand
from Bio import SeqIO

__version__ = str(Version.coerce(get_distribution('sunbeam').version))

def build_sample_list(samplelist_fp, paired_end):
    """
    Build a list of samples from a sample list file.
    :param samplelist_fp: a Path to a whitespace-delimited samplelist file,
       where the first entry is the sample name and the rest is ignored.
    :returns: A dictionary of samples with sample name and associated file(s)
    """
    Samples = {}
    with open(str(samplelist_fp)) as f:
        reader = csv.DictReader(f, fieldnames=['sample','1','2'])
        for row in reader:
            sample = row['sample']
            try:
                r1 = _verify_path(row['1'])
            except ValueError:
                raise ValueError("Associated file for {} not found.".format(
                    sample))
            r2 = None
            if paired_end:
                try:
                    r2 = _verify_path(row['2'])
                except ValueError:
                    raise ValueError(
                        "Paired-end files specified, but mate pair for '{}' "
                        "is missing or does not exist.".format(sample))
            Samples[sample] = {'1': r1, '2': r2}
    return Samples

def guess_format_string(fnames, paired_end=True, split_pattern="([_\.])"):
    """
    Try to guess the format string given a list of filenames.
    :param fnames: a list of filename strings
    :param paired_end: if True, will try to find a read-pair designator
    :param split_pattern: regex to split filenames on
    :returns: a format string with a {sample} and (maybe) an {rp} element
    """
    
    if isinstance(fnames, str):
        raise ValueError("Need a list of filenames, not a string")
    if len(fnames) == 1:
        raise ValueError("Need a list of filenames, not just one")
    if len(set(fnames)) == 1:
        raise ValueError("All filenames are the same")
    
    splits = [re.split(split_pattern, fname) for fname in fnames]

    if len(set([len(p) for p in splits])) > 1:
        raise ValueError("Files have inconsistent numbers of _ or . characters")

    elements = []
    variant_idx = []

    for i, parts in enumerate(zip(*splits)):
        items = set(parts)
        # If they're all the same, it's a common part; so add it to the element
        # list unchanged
        if len(items) == 1:
            elements.append(parts[0])
        else:
            if paired_end:
                # If all the items in a split end with 1 or 2, and only have
                # one char preceding that that's the same among all items,
                # then it's like a read-pair identifier.
                if set(_[-1] for _ in items) == {'1', '2'}:
                    prefixes = set(_[:-1] for _ in items)
                    if len(prefixes) == 1 and all(len(p) == 1 for p in prefixes):
                        prefix = parts[0][:-1]
                        elements.append(prefix)
                        elements.append("{rp}")
                        continue
            variant_idx.append(i)                        
            elements.append("{sample}")
            
    # Combine multiple variant elements
    _min = min(variant_idx)
    _max = max(variant_idx)
    elements[_min:_max+1] = ["{sample}"]
    return "".join(elements)



def _verify_path(fp):
    path = Path(fp)
    if not path.is_file():
        raise ValueError("File not found")
    return str(path.resolve())

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
