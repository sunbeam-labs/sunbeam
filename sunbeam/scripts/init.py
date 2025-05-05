import argparse
import shutil
import sys
from pathlib import Path
from sunbeam import CONFIGS_DIR
from sunbeam.project.sample_list import SampleList
from sunbeam.project.sunbeam_config import SunbeamConfig
from sunbeam.project.sunbeam_profile import SunbeamProfile


def main(argv=sys.argv):
    parser = main_parser()
    args = parser.parse_args(argv)

    project_fp = Path(args.project_fp)
    if project_fp.exists() and not args.force:
        raise SystemExit("Error: project folder already exists. Use -f to overwrite.")
    if args.force:
        shutil.rmtree(project_fp, ignore_errors=True)

    project_fp.mkdir(parents=True, exist_ok=True)

    profile_fp = CONFIGS_DIR / f"{args.profile}_profile.yaml"
    profile = SunbeamProfile.from_template(profile_fp)
    profile.to_file(project_fp / "config.yaml")

    config_fp = args.template if args.template else CONFIGS_DIR / "default_config.yml"
    config = SunbeamConfig.from_template(config_fp, project_fp)
    if args.single_end:
        config.config["all"]["paired_end"] = False
    config.to_file(project_fp / "sunbeam_config.yml")

    data_fp = args.data_fp
    if data_fp:
        data_fp = Path(data_fp)
        paired_end = not args.single_end
        sample_list = SampleList(data_fp, paired_end, args.format)
        sample_list.to_file(project_fp / "samples.csv")


def main_parser():
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

    configs = parser.add_argument_group("config file options")
    configs.add_argument(
        "--template",
        metavar="PATH",
        help="custom config file template, in YAML format",
        type=Path,
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

    return parser
