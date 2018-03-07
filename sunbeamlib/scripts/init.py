import os
import sys
import argparse
import ruamel.yaml

from pathlib import Path

from sunbeamlib import config

def _find_conda_fp():
    """Try to detect conda install in the path."""
    try:
        path = os.environ["PATH"].split(":")
        conda_fp = Path([p for p in path if "conda" in p][0]).parent
        return conda_fp
    except (KeyError, IndexError):
        pass

    
def main(argv):
    """Create a blank config file with the given name."""

    conda_fp = _find_conda_fp()
     
    parser = argparse.ArgumentParser(
        "init", description="Initialize a new sunbeam project")
    parser.add_argument("project_fp", help="Project root directory")
    parser.add_argument(
        "--defaults", choices=['microb120','pmacs','testing'],
        help="Some servers and configurations have prebuilt default configs")
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

    
    cfg = config.new(
        conda_fp=args.conda_fp, project_fp=args.project_fp, template=args.template)
    if args.defaults:
        defaults = config.load_defaults(args.defaults)
        cfg = config.update(cfg, defaults)
        
    config.dump(cfg)
