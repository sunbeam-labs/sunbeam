import os
import sys
import argparse

from pathlib import Path

from sunbeamlib.config import create_blank_config

def _find_conda_fp():
    """Try to detect conda install in the path."""
    try:
        path = os.environ["PATH"].split(":")
        conda_fp = Path([p for p in path if "conda" in p][0]).parent
        return conda_fp
    except (KeyError, IndexError):
        pass

    
def main():
    """Create a blank config file with the given name."""

    conda_fp = _find_conda_fp()
     
    parser = argparse.ArgumentParser(
        "init", description="Initialize a new sunbeam project")
    parser.add_argument("project_fp", help="Project root directory")
    parser.add_argument(
        "--server", choices=['microb120','microb191','pmacs','respublica'],
        help="Some servers have prebuilt configs; specify one here to use it.")
    parser.add_argument(
        "--conda_fp", default=conda_fp, type=Path,
        help="path to conda (default: %(default)s)"
    )
    parser.add_argument(
        "--custom_configfile_fp", help="Custom config file absolute directory", 
        type=argparse.FileType("r"), default=None)

    args = parser.parse_args()
    
    if not args.conda_fp or not args.conda_fp.exists():
        raise SystemExit((
            "Error: conda installation could not be found at '{conda_fp}'. " 
            "Specify a valid path to conda using --conda_fp."
        ).format(conda_fp = args.conda_fp if args.conda_fp else ""))
    else:
        args.conda_fp = args.conda_fp.absolute()

    if args.custom_configfile_fp is None:
        config = create_blank_config(
            args.conda_fp, args.project_fp, template=args.server)
    elif os.path.isfile(args.custom_configfile_fp.name):
        with open(args.custom_configfile_fp.name, 'r') as f:
            config = f.read()
    else:
        raise SystemExit(
            "Specify a valid path to custom config file --custom_configfile_fp.")
 
    sys.stdout.write(config)

