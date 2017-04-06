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
        "--template", default=None, 
        help="Path to custom config file template, in YAML format", 
        type=argparse.FileType("r"))

    args = parser.parse_args()
    
    if not args.conda_fp or not args.conda_fp.exists():
        raise SystemExit((
            "Error: conda installation could not be found at '{conda_fp}'. " 
            "Specify a valid path to conda using --conda_fp."
        ).format(conda_fp = args.conda_fp if args.conda_fp else ""))
    else:
        args.conda_fp = args.conda_fp.absolute()

    if args.template is None:
        config = create_blank_config(
            args.conda_fp, args.project_fp, template=args.server)
    else:
        config = args.template.read().format(
                 CONDA_FP=args.conda_fp, PROJECT_FP=args.project_fp)
 
    sys.stdout.write(config)
