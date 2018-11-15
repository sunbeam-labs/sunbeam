import os
import sys
import argparse
import ruamel.yaml
from pathlib import Path

from sunbeamlib import config
from .init import check_existing
import subprocess
import io
import re
    
def main(argv=sys.argv):

    try:
        conda_prefix = os.environ.get("CONDA_PREFIX")
    except (KeyError, IndexError):
        raise SystemExit(
            "Could not determine Conda prefix. Activate your Sunbeam "
            "environment and try this command again.")

    description_str = (
        "Initialize a new Sunbeam project in a given directory, creating "
        "a new config file and (optionally) a sample list.")
    
    parser = argparse.ArgumentParser(
        "init", description=description_str)
    parser.add_argument(
        "project_fp", type=Path,
        help="project directory (will be created if it does not exist)")
    parser.add_argument(
        "-f", "--force", help="overwrite files if they already exist",
        action="store_true")
    parser.add_argument(
        "--output", help=(
            "name of config file (%(default)s)"),
        default="sunbeam_config.yml", metavar="FILE")

    configs = parser.add_argument_group("config file options")
    configs.add_argument(
        "--defaults", type=argparse.FileType("r"), metavar="FILE",
        help="set of default values to use to populate config file")
    configs.add_argument(
        "--template", default=None, metavar="FILE",
        help="custom config file template, in YAML format", 
        type=argparse.FileType("r"))

    samplelist = parser.add_argument_group("sample list options")
    samplelist.add_argument(
        "--sra_acc", metavar="PATH",
        help="list of SRA accession numbers")

    args = parser.parse_args(argv)

    # Create project folder if it doesn't exist
    project_fp_exists = False
    project_fp = args.project_fp
    
    try:
        project_fp = args.project_fp.resolve()
        project_fp_exists = project_fp.exists()
    except FileNotFoundError:
        pass
    
    if not project_fp_exists:
        sys.stderr.write(
            "Creating project folder at {}...\n".format(args.project_fp))
        project_fp.mkdir(parents=True, exist_ok=True)
    

    # Check if files already exist
    config_file = check_existing(project_fp/args.output, args.force)
    samplelist_file = check_existing(project_fp/"samples.csv", args.force)

    # Create config file        
    cfg = config.new(
        conda_fp=conda_prefix,
        project_fp=project_fp,
        template=args.template)

    if args.defaults:
        defaults = ruamel.yaml.safe_load(args.defaults)
        cfg = config.update(cfg, defaults)

    with config_file.open('w') as out:
        config.dump(cfg, out)

    sys.stderr.write("New config file written to {}\n".format(config_file))

    # Create sample list
    try:
        with samplelist_file.open('w') as out:
            build_sample_list_sra(accessions = args.data_acc output_file = out)
            sys.stderr.write(
                "New sample list written to {}\n".format(samplelist_file))

def build_sample_list_sra(accessions, output_file):
    """
    Create samples CSV file with all samples from given SRA accessions.

    :param accessions: list of SRA accession strings
    :param output_file: path to CSV file to write
    """

    samples = find_samples_sra(accessions)

    # How many rows do we have for each entry?  Simple case is just one or two,
    # but it's possible we'll have both.
    lengths = {len(fqs) for fqs in samples.values()}
    fmt = "Found {} samples, {}.\n"
    if lengths == {1}:
        # OK, unpaired.
        is_single_end = True
        sys.stderr.write(fmt.format(len(samples), ["paired", "unpaired"][is_single_end]))
        # TODO write CSV
    elif lengths == {2}:
        # OK, paired-end.
        is_single_end = False
        sys.stderr.write(fmt.format(len(samples), ["paired", "unpaired"][is_single_end]))
        # TODO write CSV
    elif lengths == {1,2}:
        # A mix of both.  Sunbeam just does one or the other so the user must
        # choose.
        # TODO raise Warning and write the two sets separately.
        sys.stderr.write("Found {} samples.\n".format(len(samples)))
        pass
    else:
        # ??? bail.
        bad_lens = [x for x in lengths if x not in [1,2]]
        bad_lens = ', '.join([str(x) for x in bad_lens])
        raise SystemExit("SRA metadata must specify 1 or 2 files per sample;"
                "not %s" % bad_lens)

    fieldnames = ["sample", "1", "2"]
    writer = csv.DictWriter(output_file, fieldnames=fieldnames)
    for sample in samples.keys():
        writer.writerow({'sample':sample, **samples[sample]})


def find_samples_sra(accessions, dir_fp="download"):
    """
    Build a dict of samples from SRA accession numbers.
    
    Accession numbers will be automatically expanded to sample-level values
    (for example, all samples belonging to a given BioProject).
    
    :param accessions: list of SRA accession strings
    :param dir_fp: directory to prepend to data file paths.
    :returns: A dictionary of samples, with sample names as keys
    """
    dir_fp = Path(dir_fp)
    # Call grabseqs to get raw text data listing accessions
    args = ["grabseqs", "sra", "-l"] + accessions
    result = subprocess.run(args, stdout=subprocess.PIPE, check=True)
    file_like = io.StringIO(result.stdout.decode("ASCII"))
    # Parse as CSV into list of lists
    reader = csv.reader(file_like)
    data = list(reader)
    # Figure out accession number for each entry (since we know what the format
    # is) and structure as a dict, with accessions as keys and file paths
    # (using given path) as values.
    samp_parse = lambda txt: re.match('^([0-9A-Za-z]+)[_\.]', txt).group(1)
    prepend_fp = lambda files: [str(dir_fp/fp) for fp in files]
    try:
        Samples = {samp_parse(d[0]): prepend_fp(d) for d in data}
    except AttributeError:
        raise SystemExit("SRA file name did not match expected pattern")
    return(Samples)
