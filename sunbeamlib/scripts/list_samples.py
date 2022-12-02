import sys
from pathlib import Path
import argparse
import csv

from snakemake.utils import listfiles

from sunbeamlib import guess_format_string, SampleFormatError, MissingMatePairError


def main(argv=sys.argv):

    parser = argparse.ArgumentParser("sunbeam list_samples")

    parser.add_argument("data_fp", help="Path to folder containing reads", type=Path)

    parser.add_argument(
        "-s",
        "--single_end",
        action="store_true",
        help="Samples are single-end (not paired-end)",
    )

    parser.add_argument(
        "-f",
        "--format",
        metavar="STR",
        help=("Filename format (e.g. {sample}_R{rp}.fastq.gz). (default: " "guessed)"),
    )

    args = parser.parse_args(argv)
    try:
        build_sample_list(args.data_fp, args.format, sys.stdout, args.single_end)
    except SampleFormatError as e:
        raise SystemExit(
            "Error: could not build sample list. Specify correct sample filename"
            " format using --format. \n\tReason: {}".format(e)
        )
    except MissingMatePairError as e:
        raise SystemExit(
            "Assuming paired-end reads, but could not find mates. Specify "
            "--single_end if not paired-end, or provide sample name format "
            "using --format."
            "\n\tReason: {}".format(e)
        )


def build_sample_list(data_fp, format_str, output_file, is_single_end):

    data_fp = data_fp.resolve()
    fnames = [f.name for f in data_fp.iterdir() if f.is_file()]

    if not format_str:
        sys.stderr.write(
            "Guessing sample name format from files in {}...\n".format(data_fp)
        )
        format_str = guess_format_string(fnames, not is_single_end)
        sys.stderr.write("  Best guess: {}\n".format(format_str))

    samples = find_samples(data_fp, format_str)

    # Check for mate pairs if single end
    if not is_single_end:
        no_match = []
        for sample, reads in samples.items():
            if "2" not in reads:
                no_match.append(sample)
        if len(no_match) > 0:
            raise MissingMatePairError(
                "missing mate pairs for samples: {} ".format(", ".join(no_match))
            )

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
    files = list(listfiles(str(data_fp / filename_fmt)))
    Samples = {t[1]["sample"]: {} for t in files}
    for f in files:
        fpath = f[0]
        wcards = f[1]
        rp = wcards.get("rp")
        if rp:
            if not rp in ["1", "2"]:
                raise ValueError(
                    "'{rp}' should capture just '1' or '2' in filename, nothing else"
                )
            Samples[wcards["sample"]][rp] = fpath
        else:
            Samples[wcards["sample"]]["1"] = fpath
    return Samples


def _write_samples_csv(samples, out):
    fieldnames = ["sample", "1", "2"]
    writer = csv.DictWriter(out, fieldnames=fieldnames)
    for sample in samples.keys():
        writer.writerow({"sample": sample, **samples[sample]})


def check_existing(path, force=False):
    if path.is_dir():
        raise SystemExit(
            "Error: specified file '{}' exists and is a directory".format(path)
        )
    if path.is_file() and not force:
        raise SystemExit(
            "Error: specified file '{}' exists. Use --force to "
            "overwrite.".format(path)
        )
    return path
