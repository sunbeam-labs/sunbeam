import os
import sys
import argparse
import ruamel.yaml
from pathlib import Path

from .list_samples import (
        build_sample_list,
        build_sample_list_sra,
        MissingMatePairError,
        SampleFormatError,
        check_existing)
from sunbeamlib import config

def main(argv=sys.argv):
    """Create a blank config file with the given name."""

    conda_prefix = get_conda_prefix()
    args = parse_args(argv)
    project_fp = setup_project_folder(args)
    samplelists = write_samples_from_input(args, project_fp)
    write_config(args, conda_prefix, project_fp, samplelists)

def get_conda_prefix():
    try:
        conda_prefix = os.environ.get("CONDA_PREFIX")
    except (KeyError, IndexError):
        raise SystemExit(
            "Could not determine Conda prefix. Activate your Sunbeam "
            "environment and try this command again.")
    return conda_prefix

def parse_args(argv):
    description_str = (
        "Initialize a new Sunbeam project in a given directory, creating "
        "a new config file and (optionally) a sample list.  The sample list "
        "source can be either a folder of input files or a list of SRA accession numbers.")
    
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

    samplelist = parser.add_argument_group("sample list options",
            ("Options to automatically generate a sample list. --data_fp (for "
             "reading filenames from a specified folder) and --data_acc (for "
             "fetching SRA accession numbers to be downloaded during a run) "
             "are mutually exclusive.  If --data_acc is used, all SRA Run "
             "(sample) entries corresponding to the given accession numbers "
             "(e.g. PRJNA###, SAMN###, SRR###) will be used to create the "
             "sample list.  Note in this case project_fp should either be "
             "given before --data_acc or separated by a '--' to distinguish "
             "it from an accession number."))
    
    samplelist.add_argument(
        "--data_fp", type=Path, metavar="PATH",
        help="path to folder containing fastq.gz files")
    samplelist.add_argument(
        "--format", type=str, metavar="STR",
        help="filename format for --data_fp (default: guessed)")
    samplelist.add_argument(
        "--single_end", action="store_true",
        help="fastq files are in single-end, not paired-end, format for --data_fp")
    samplelist.add_argument("--data_acc", metavar="ACC", nargs="+",
        help="list of SRA-compatible accession numbers")

    # argparse doesn't support complete argument groups that are mutually
    # exclusive (need to do subparsers/subcommands) but this seems good enough:
    # https://stackoverflow.com/a/27675614
    args = parser.parse_args(argv)
    if args.data_fp and args.data_acc:
        parser.error("--data_fp and --data_acc are mutually exclusive options")
    return args

def setup_project_folder(args):
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
    return project_fp

def write_samples_from_input(args, project_fp):
    """Write sample list CSV from existing files."""
    if args.data_acc:
        # SRA case: create one (for all unpaired or paired) or two (if both
        # present) samples CSV files.
        samplelists = build_sample_list_sra(
                accessions = args.data_acc,
                project_fp = project_fp,
                force = args.force)
        for fp in samplelists.values():
            sys.stderr.write("New sample list written to {}\n".format(fp))
    else:
        # filnames from local disk case: create one samples CSV file.
        samplelist_file = check_existing(project_fp/"samples.csv", args.force)
        if args.data_fp:
            try:
                with samplelist_file.open("w") as out:
                    build_sample_list(
                        data_fp = args.data_fp,
                        format_str = args.format,
                        output_file = out,
                        is_single_end = args.single_end)
            except SampleFormatError as e:
                raise SystemExit(
                    "Error: could not create sample list. Specify correct sample filename"
                    " format using --format.\n  Reason: {}".format(e))
            except MissingMatePairError as e:
                raise SystemExit(
                    "Error: assuming paired-end reads, but could not find mates. Specify "
                    "--single-end if not paired-end, or provide sample name format "
                    "using --format."
                    "\n  Reason: {}".format(e))
            sys.stderr.write(
                "New sample list written to {}\n".format(samplelist_file))
        samplelists = {["paired", "unpaired"][args.single_end]: samplelist_file}
    return samplelists

def write_config(args, conda_prefix, project_fp, samplelists):
    multiple_configs = len(samplelists.values()) != 1
    for layout in samplelists.keys(): # one paired, one unpaired, or one of each
        if multiple_configs:
            config_file = check_existing(project_fp/Path(layout+"_"+str(Path(project_fp/args.output).name)), args.force)
        else:
            config_file = check_existing(project_fp/args.output, args.force)
        cfg = config.new(
            conda_fp=conda_prefix,
            project_fp=project_fp,
            template=args.template)
        defaults = {}
        if args.defaults:
            defaults = ruamel.yaml.safe_load(args.defaults)
        # Override loaded config defaults (if any) for a few specific items.
        paired = layout == "paired"
        defaults["all"] = defaults.get("all", {})
        defaults["all"]["paired_end"] = paired
        defaults["all"]["samplelist_fp"] = samplelists[layout].name
        if args.data_acc:
            defaults["all"]["download_reads"] = True
        # Convert cfg from raw text to dict, including any defaults we have
        # set, and write to disk.
        cfg = config.update(cfg, defaults)
        with config_file.open('w') as out:
            config.dump(cfg, out)
        sys.stderr.write("New config file written to {}\n".format(config_file))
    if multiple_configs:
        raise SystemExit("Found both paired and unpaired reads. Wrote two sample lists "
                        "and config files, with '_paired' or '_single' appended.")
