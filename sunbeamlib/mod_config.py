import sys
import ruamel.yaml
import argparse

from .config import update_config

def main():

    parser = argparse.ArgumentParser(
        description="Modifies a Sunbeam config file with either a YAML file or string")
    parser.add_argument(
        "--config", help="Config file to modify", type=argparse.FileType("r"),
        default=sys.stdin)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--mod_fp", help="YAML file with updated config values",
        type=argparse.FileType("r"), metavar="FILE")
    group.add_argument(
        "--mod_str", help="YAML string (e.g. 'blast:{threads:4}')",
        type=str, metavar="STR")

    args = parser.parse_args()

    if args.mod_fp:
        mods = ruamel.yaml.safe_load(args.mod_fp)
    else:
        mods = ruamel.yaml.safe_load(args.mod_str)
    print(mods)
    config = update_config(args.config, mods)
    ruamel.yaml.round_trip_dump(config, sys.stdout)

if __name__ == "__main__":
    main()
