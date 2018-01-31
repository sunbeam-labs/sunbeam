import sys
import ruamel.yaml
import argparse

from ..config import update

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Modifies a Sunbeam config file with a YAML file or string, or "
            "values from an old config"))
    
    parser.add_argument(
        "--config", help="Config file to modify", type=argparse.FileType("r"),
        default=sys.stdin)

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument(
        "--old_config", help="Old config file",
        type=argparse.FileType("r"), metavar="FILE")
    group.add_argument(
        "--mod_fp", help="YAML file with updated config values",
        type=argparse.FileType("r"), metavar="FILE")
    group.add_argument(
        "--mod_str", help="YAML string (e.g. 'blast: {threads: 4}')",
        type=str, metavar="STR")

    args = parser.parse_args()

    # Use "strict" updating if sourcing from old config (i.e. do not make any
    # new keys)
    strict = False
    if args.mod_fp:
        mods = ruamel.yaml.safe_load(args.mod_fp)
    elif args.old_config:
        strict = True
        mods = ruamel.yaml.safe_load(args.old_config)
        # Remove version number from old config so we don't transfer it over
        mods.get('all', {}).pop('version', None)
    else:
        mods = ruamel.yaml.safe_load(args.mod_str)
        if isinstance(mods, str):
            sys.exit(
                "Invalid YAML: did you make sure to put spaces between keys and "
                "values?")

    config = update(args.config, mods, strict)
    ruamel.yaml.round_trip_dump(config, sys.stdout)


    
