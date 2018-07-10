import os
import sys
import argparse
import ruamel.yaml
from pathlib import Path

from .list_samples import build_sample_list, MissingMatePairError, SampleFormatError
from sunbeamlib import config
    
def main(argv=sys.argv):
    """Create a blank config file with the given name."""

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
        "--data_fp", type=Path, metavar="PATH",
        help="path to folder containing fastq.gz files")
    samplelist.add_argument(
        "--format", type=str, metavar="STR",
        help="filename format (default: guessed)")
    samplelist.add_argument(
        "--single_end", action="store_true",
        help="fastq files are in single-end, not paired-end, format")

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

    # Create sample list (if given data_fp)
    if args.data_fp:
        try:
            with samplelist_file.open('w') as out:
                build_sample_list(
                    data_fp = args.data_fp,
                    format_str = args.format,
                    output_file = out,
                    is_single_end = args.single_end)
                sys.stderr.write(
                    "New sample list written to {}\n".format(samplelist_file))
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
        
def check_existing(path, force=False):
    if path.is_dir():
        raise SystemExit(
            "Error: specified file '{}' exists and is a directory".format(path))
    if path.is_file() and not force:
        raise SystemExit(
            "Error: specified file '{}' exists. Use --force to "
            "overwrite.".format(path))
    return path
