import sys
from pathlib import Path
import argparse
import csv

import ruamel.yaml
from snakemake.utils import listfiles

from sunbeamlib import guess_format_string

def main(argv=sys.argv):

    parser = argparse.ArgumentParser(
        "sunbeam list_samples"
    )

    parser.add_argument(
        "data_fp", help="Path to folder containing reads", type=Path)

    parser.add_argument(
        "-s", "--single_end", action="store_true",
        help="Samples are single-end (not paired-end)")

    parser.add_argument(
        "-f", "--format", metavar="STR", help=(
            "Filename format (e.g. {sample}_R{rp}.fastq.gz). (default: "
            "guessed)"))

    args = parser.parse_args(argv)

    data_fp = args.data_fp.resolve()
    fnames = [f.name for f in data_fp.iterdir() if f.is_file()]
    
    if not args.format:
        try:
            sys.stderr.write("Attempting to guess sample name format... ")
            args.format = guess_format_string(fnames, not args.single_end)
            sys.stderr.write("'{}'\n".format(args.format))
        except ValueError as e:
            raise SystemExit(
                "Could not guess sample name format. Please provide it using"
                " --format. (Reason: {})".format(e))

    samples = find_samples(args.data_fp.resolve(), args.format)

    if len(samples) == 0:
        raise SystemExit("Error: No samples matching the given format found.")

    sys.stderr.write("Found {} samples matching this format.\n".format(len(samples)))
    fieldnames = ["sample", "1", "2"]
    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
    for sample in samples.keys():
        writer.writerow({'sample':sample, **samples[sample]})

def find_samples(data_fp, filename_fmt):
    """
    Build a list of samples from a data filepath and filename format.

    :param data_fp: a Path to data files
    :param filename_fmt: a string giving wildcards for {sample} and (optionally)
        {rp}, for read pair (e.g. R1 or R2).
    :returns: A dictionary of samples, with sample names as keys:
       Samples = {
         'sample1': {
           '1': 'path/to/sample1_R1.fastq.gz',
           '2': 'path/to/sample1_R2.fastq.gz',
           'paired_end': True
         }, ...
       }
    """
    files = list(listfiles(str(data_fp/filename_fmt)))
    Samples = {t[1]['sample']: {} for t in files}
    for f in files:
        fpath = f[0]
        wcards = f[1]
        rp = wcards.get('rp')
        if rp:
            if not rp in ['1', '2']:
                raise ValueError(
                    "'{rp}' should capture just '1' or '2' in filename, nothing else")
            Samples[wcards['sample']][rp] = fpath
        else:
            Samples[wcards['sample']]['1'] = fpath
    return Samples

        
