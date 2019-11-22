import sys
import io
import re
from pathlib import Path
import argparse
import csv
import subprocess

import ruamel.yaml
from snakemake.utils import listfiles

from sunbeamlib import guess_format_string, SampleFormatError, MissingMatePairError

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
    try:
        build_sample_list(args.data_fp, args.format, sys.stdout, args.single_end)
    except SampleFormatError as e:
        raise SystemExit(
            "Error: could not build sample list. Specify correct sample filename"
            " format using --format. \n\tReason: {}".format(e))
    except MissingMatePairError as e:
        raise SystemExit(
            "Assuming paired-end reads, but could not find mates. Specify "
            "--single_end if not paired-end, or provide sample name format "
            "using --format."
            "\n\tReason: {}".format(e))

def build_sample_list(data_fp, format_str, output_file, is_single_end):

    data_fp = data_fp.resolve()
    fnames = [f.name for f in data_fp.iterdir() if f.is_file()]
    
    if not format_str:
        sys.stderr.write(
            "Guessing sample name format from files in {}...\n".format(data_fp))
        format_str = guess_format_string(fnames, not is_single_end)
        sys.stderr.write("  Best guess: {}\n".format(format_str))
        
    samples = find_samples(data_fp, format_str)

    # Check for mate pairs if single end
    if not is_single_end:
        no_match = []
        for sample, reads in samples.items():
            if '2' not in reads:
                no_match.append(sample)
        if len(no_match) > 0:
            raise MissingMatePairError(
                "missing mate pairs for samples: {} ".format(
                    ", ".join(no_match)))

    if len(samples) == 0:
        raise SampleFormatError("no samples matching the given format found.")

    sys.stderr.write("Found {} samples in {}.\n".format(len(samples), data_fp))
    _write_samples_csv(samples, output_file)

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
           '2': 'path/to/sample1_R2.fastq.gz', #optional
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

def build_sample_list_sra(accessions, project_fp, force):
    """
    Create samples CSV file with all samples from given SRA accessions.

    :param accessions: list of SRA accession strings
    :param project_fp: project folder path
    :param force: force file creation?
    """

    samples = find_samples_sra(accessions)

    # How many rows do we have for each entry?  Simple case is just one or two,
    # but it's possible we'll have both.
    lengths = {len(fqs) for fqs in samples.values()}
    fmt = "Found {} samples, {}.\n"
    cases = ["unpaired", "paired"]
    files = {}
    if lengths == {1} or lengths == {2}:
        # OK, either all paired or all unpaired.
        case = cases[list(lengths)[0] - 1]
        sys.stderr.write(fmt.format(len(samples), case))
        fp = check_existing(project_fp/"samples.csv", force)
        files[case] = fp
        with fp.open("w") as f:
            _write_samples_csv(samples, f)
    elif lengths == {1,2}:
        # A mix of both.  Sunbeam just does one or the other so the user must
        # choose.
        sys.stderr.write(fmt.format(len(samples), "both paired and unpaired"))
        sys.stderr.write("config file and samples will be available with both"
                " _paired and _unpaired suffixes.  These can be run separately"
                " with sunbeam run by passing the appropriate config file.\n")
        for case in cases:
            fp = project_fp/("samples_%s.csv" % case)
            fp = check_existing(fp, force)
            len_exp = cases.index(case) + 1
            samps = {k: v for k, v in samples.items() if len(v) == len_exp}
            with fp.open("w") as f:
                _write_samples_csv(samps, f)
            files[case] = fp
    else:
        # ??? bail.
        bad_lens = [x for x in lengths if x not in [1,2]]
        bad_lens = ', '.join([str(x) for x in bad_lens])
        raise SystemExit("SRA metadata must specify 1 or 2 files per sample;"
                "not %s" % bad_lens)
    return(files)

def find_samples_sra(accessions, dir_fp="download"):
    """
    Build a dict of samples from SRA accession numbers.
    
    Accession numbers will be automatically expanded to sample-level values
    (for example, all samples belonging to a given BioProject).
    The returned dict will be structured the same way as those created by
    find_samples above.
    
    :param accessions: list of SRA accession strings
    :param dir_fp: directory to prepend to data file paths.
    :returns: A dictionary of samples, with sample names as keys
    """
    dir_fp = Path(dir_fp)
    # Call grabseqs to get raw text data listing accessions
    cmd_args = ["grabseqs", "sra", "-l"] + accessions
    result = subprocess.run(cmd_args, stdout=subprocess.PIPE, check=True)
    file_like = io.StringIO(result.stdout.decode("ASCII"))
    # Parse as CSV into list of lists
    reader = csv.reader(file_like)
    data = list(reader)
    # Figure out accession number for each entry (since we know what the format
    # is) and structure as a dict, with accessions as keys and file paths
    # (using given path) as values.
    samp_parse = lambda txt: re.match('^([0-9A-Za-z]+)[_\.]', txt).group(1)
    prepend_fp = lambda files: [str(dir_fp/fp) for fp in files]
    dictify = lambda files: {str(i+1): v for i, v in zip(range(len(files)), files)}
    try:
        Samples = {samp_parse(d[0]): dictify(prepend_fp(d)) for d in data}
    except AttributeError:
        raise SystemExit("SRA file name did not match expected pattern")
    return(Samples)

def _write_samples_csv(samples, out):
    fieldnames = ["sample", "1", "2"]
    writer = csv.DictWriter(out, fieldnames=fieldnames)
    for sample in samples.keys():
        writer.writerow({'sample':sample, **samples[sample]})

def check_existing(path, force=False):
    if path.is_dir():
        raise SystemExit(
            "Error: specified file '{}' exists and is a directory".format(path))
    if path.is_file() and not force:
        raise SystemExit(
            "Error: specified file '{}' exists. Use --force to "
            "overwrite.".format(path))
    return path
