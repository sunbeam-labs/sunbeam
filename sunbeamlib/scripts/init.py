import os
import sys
import argparse
import ruamel.yaml
from pathlib import Path

from sunbeamlib import config
    
def main(argv=sys.argv):
    """Create a blank config file with the given name."""

    try:
        conda_prefix = os.environ.get("CONDA_PREFIX")
    except (KeyError, IndexError):
        raise SystemExit(
            "Could not determine Conda prefix. Activate your Sunbeam "
            "environment and try this command again.")
     
    parser = argparse.ArgumentParser(
        "init", description="Initialize a new sunbeam project")
    parser.add_argument(
        "-o", "--output", help="Name of config file (written to project_fp)",
        default="sunbeam_config.yml", metavar="FILE")
    parser.add_argument(
        "-f", "--force", help="Overwrite config file if it already exists",
        action="store_true")
    parser.add_argument(
        "-d", "--defaults", type=argparse.FileType("r"), metavar="FILE",
        help="Set of default values to use to populate config file")
    parser.add_argument(
        "project_fp", type=Path,
        help="Project directory (will be created if it does not exist)")
    parser.add_argument(
        "--template", default=None, metavar="FILE",
        help="Custom config file template, in YAML format", 
        type=argparse.FileType("r"))

    args = parser.parse_args(argv)

    # Create project folder if it doesn't exist
    if not args.project_fp.exists():
        sys.stderr.write("Creating project folder at {}...\n".format(args.project_fp))
        args.project_fp.mkdir(parents=True, exist_ok=True)

    cfg_fp = args.project_fp/args.output
    if cfg_fp.exists() and not args.force:
        raise SystemExit(
            "Error: config file already exists at {}. Choose a new name or use --force "
            "to overwrite.".format(cfg_fp))
        
    cfg = config.new(
        conda_fp=conda_prefix,
        project_fp=args.project_fp,
        template=args.template)

    if args.defaults:
        defaults = ruamel.yaml.safe_load(args.defaults)
        cfg = config.update(cfg, defaults)

    with open(cfg_fp, 'w') as out:
        config.dump(cfg, out)

    sys.stderr.write("New config file written to {}\n".format(cfg_fp))
