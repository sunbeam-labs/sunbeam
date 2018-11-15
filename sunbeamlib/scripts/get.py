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
import csv
    
def main(argv=sys.argv):

    try:
        conda_prefix = os.environ.get("CONDA_PREFIX")
    except (KeyError, IndexError):
        raise SystemExit(
            "Could not determine Conda prefix. Activate your Sunbeam "
            "environment and try this command again.")

    description_str = (
        "Initialize a new Sunbeam project in a given directory, creating "
        "a new config file and sample list from SRA accessions.")

    parser = argparse.ArgumentParser(
        "get", description=description_str)
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
    samplelist.add_argument("--data_acc", metavar="PATH", nargs="+",
        help="list of SRA accession numbers")

    args = parser.parse_args(argv)

    ### Project Directory

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

    ### Sample list(s)

    # Create sample list(s).  One or two depending on paired/unpaired
    # situation.
    samplelists = build_sample_list_sra(args.data_acc, args)
    for fp in samplelists.values():
        sys.stderr.write("New sample list written to {}\n".format(fp))

    ### Config(s)

    # TODO from here on, update config-handling to handle possible mixed
    # unpaired/paired case.  (Append suffix to existing args.output in that
    # case, between name and .yml?)

    # Louis thought: gonna append for now since I don't know how to parse Path objects

    multiple_configs = len(samplelists.values()) != 1
    for layout in samplelists.keys(): # one paired, one unpaired, or one of each
        # Create config file
        if multiple_configs:
            config_file = check_existing(str(project_fp/args.output)+"_"+layout, args.force)
        else:
            config_file = check_existing(project_fp/args.output, args.force)
        cfg = config.new(
            conda_fp=conda_prefix,
            project_fp=project_fp,
            template=args.template)
        if args.defaults:
            defaults = ruamel.yaml.safe_load(args.defaults)
            cfg = config.update(cfg, defaults)
        if layout == "paired":
            paired = "true"
        else:
            paired = "false"
        sample_file = samplelists[layout].name
        cfg = config.update(cfg, {"all":{"paired":paired, "download_reads":"true"}})
        with config_file.open('w') as out:
            config.dump(cfg, out)
        sys.stderr.write("New config file written to {}\n".format(config_file))
    if multiple_configs:
        raise Exception("Found both paired and unpaired reads. Wrote two sample lists "
                        "and config files, with '_paired' or '_single' appended.")

def _write_samples_csv(samples, fp):
    fieldnames = ["sample", "1", "2"]
    with fp.open('w') as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        for sample in samples.keys():
            writer.writerow({'sample':sample, **samples[sample]})

def build_sample_list_sra(accessions, args):
    """
    Create samples CSV file with all samples from given SRA accessions.

    :param accessions: list of SRA accession strings
    :param args: arguments from command-line
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
        fp = check_existing(args.project_fp/"samples.csv", args.force)
        files[case] = fp
        _write_samples_csv(samples, fp)
    elif lengths == {1,2}:
        # A mix of both.  Sunbeam just does one or the other so the user must
        # choose.
        sys.stderr.write(fmt.format(len(samples), "both paired and unpaired"))
        sys.stderr.write("config file and samples will be available with both"
                " _paired and _unpaired suffixes.  These can be run separately"
                " with sunbeam run by passing the appropriate config file.\n")
        for case in cases:
            fp = args.project_fp/("samples_%s.csv" % case)
            fp = check_existing(fp, args.force)
            len_exp = cases.index(case) + 1
            samps = {k: v for k, v in samples.items() if len(v) == len_exp}
            _write_samples_csv(samps, fp)
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
    dictify = lambda files: {str(i+1): v for i, v in zip(range(len(files)), files)}
    try:
        Samples = {samp_parse(d[0]): dictify(prepend_fp(d)) for d in data}
    except AttributeError:
        raise SystemExit("SRA file name did not match expected pattern")
    return(Samples)
