import os
import shutil
import sys
import argparse
import ruamel.yaml
from pathlib import Path

from .list_samples import (
    build_sample_list,
    MissingMatePairError,
    SampleFormatError,
    check_existing,
)
from sunbeamlib import config


def main(argv=sys.argv):
    """Create a blank config file with the given name."""

    args = parse_args(argv)
    project_fp = setup_project_folder(args)
    samplelists = write_samples_from_input(args, project_fp)
    write_config(args, project_fp, samplelists)
    write_profile(args, project_fp)


def get_conda_prefix():
    try:
        conda_prefix = os.environ.get("CONDA_PREFIX")
    except (KeyError, IndexError):
        raise SystemExit(
            "Could not determine Conda prefix. Activate your Sunbeam "
            "environment and try this command again."
        )
    return conda_prefix


def parse_args(argv):
    description_str = (
        "Initialize a new Sunbeam project in a given directory, creating "
        "a new config file and (optionally) a sample list."
    )

    parser = argparse.ArgumentParser("init", description=description_str)
    parser.add_argument(
        "project_fp",
        type=Path,
        help="project directory (will be created if it does not exist)",
    )
    parser.add_argument(
        "-f",
        "--force",
        help="overwrite files if they already exist",
        action="store_true",
    )
    parser.add_argument(
        "--output",
        help=("name of config file (%(default)s)"),
        default="sunbeam_config.yml",
        metavar="FILE",
    )

    configs = parser.add_argument_group("config file options")
    configs.add_argument(
        "--defaults",
        type=argparse.FileType("r"),
        metavar="FILE",
        help="set of default values to use to populate config file",
    )
    configs.add_argument(
        "--template",
        default=None,
        metavar="FILE",
        help="custom config file template, in YAML format",
        type=argparse.FileType("r"),
    )

    samplelist = parser.add_argument_group("sample list options")
    samplelist.add_argument(
        "--data_fp",
        type=Path,
        metavar="PATH",
        help="path to folder containing fastq.gz files",
    )
    samplelist.add_argument(
        "--format",
        type=str,
        metavar="STR",
        help="filename format for --data_fp (default: guessed)",
    )
    samplelist.add_argument(
        "--single_end",
        action="store_true",
        help="fastq files are in single-end, not paired-end, format for --data_fp",
    )

    profile = parser.add_argument_group("profile file options")
    profile.add_argument(
        "--profile",
        type=str,
        metavar="STR",
        default="default",
        help="name of the profile template to use (e.g. default, slurm) (default: default)",
    )

    # argparse doesn't support complete argument groups that are mutually
    # exclusive (need to do subparsers/subcommands) but this seems good enough:
    # https://stackoverflow.com/a/27675614
    args = parser.parse_args(argv)
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
        sys.stderr.write("Creating project folder at {}...\n".format(args.project_fp))
        project_fp.mkdir(parents=True, exist_ok=True)
    return project_fp


def write_samples_from_input(args, project_fp):
    """Write sample list CSV from existing files."""
    samplelist_file = check_existing(project_fp / "samples.csv", args.force)
    if args.data_fp:
        try:
            with samplelist_file.open("w") as out:
                build_sample_list(
                    data_fp=args.data_fp,
                    format_str=args.format,
                    output_file=out,
                    is_single_end=args.single_end,
                )
        except SampleFormatError as e:
            os.remove(project_fp / "samples.csv")
            raise SystemExit(
                "Error: could not create sample list. Specify correct sample filename"
                " format using --format.\n  Reason: {}".format(e)
            )
        except MissingMatePairError as e:
            os.remove(project_fp / "samples.csv")
            raise SystemExit(
                "Error: assuming paired-end reads, but could not find mates. Specify "
                "--single_end if not paired-end, or provide sample name format "
                "using --format."
                "\n  Reason: {}".format(e)
            )
        sys.stderr.write("New sample list written to {}\n".format(samplelist_file))
    samplelists = {["paired", "unpaired"][args.single_end]: samplelist_file}
    return samplelists


def write_config(args, project_fp, samplelists):
    multiple_configs = len(samplelists.values()) != 1
    for layout in samplelists.keys():  # one paired, one unpaired, or one of each
        if multiple_configs:
            config_file = check_existing(
                project_fp
                / Path(layout + "_" + str(Path(project_fp / args.output).name)),
                args.force,
            )
        else:
            config_file = check_existing(project_fp / args.output, args.force)
        cfg = config.new(project_fp=project_fp, template=args.template)
        defaults = {}
        if args.defaults:
            defaults = ruamel.yaml.safe_load(args.defaults)
        # Override loaded config defaults (if any) for a few specific items.
        paired = layout == "paired"
        defaults["all"] = defaults.get("all", {})
        defaults["all"]["paired_end"] = paired
        defaults["all"]["samplelist_fp"] = samplelists[layout].name
        cfg = config.update(cfg, defaults)
        with config_file.open("w") as out:
            config.dump(cfg, out)
        sys.stderr.write("New config file written to {}\n".format(config_file))
    if multiple_configs:
        raise SystemExit(
            "Found both paired and unpaired reads. Wrote two sample lists "
            "and config files, with '_paired' or '_single' appended."
        )


def write_profile(args, project_fp):
    sunbeam_dir = Path(os.getenv("SUNBEAM_DIR", os.getcwd()))
    template_fp = sunbeam_dir / "sunbeamlib" / "data" / f"{args.profile}_profile.yaml"
    config_fp = project_fp / "config.yaml"
    shutil.copyfile(template_fp, config_fp)
    with open(config_fp, "a") as f:
        f.write("\n# Filepath of this projects configfile")
        f.write(f"\nconfigfile: {project_fp/args.output}")
    sys.stderr.write("New profile file written to {}\n".format(config_fp))
