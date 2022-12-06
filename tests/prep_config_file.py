import sys
import ruamel.yaml
import argparse
import collections


def update_dict(target, new):
    for k, v in new.items():
        if isinstance(v, collections.Mapping):
            target[k] = update_dict(target.get(k, {}), v)
        else:
            target[k] = v
    return target


def main():

    parser = argparse.ArgumentParser(description="Modifies a Sunbeam config file")
    parser.add_argument(
        "--config",
        help="Config file to modify",
        type=argparse.FileType("r"),
        default=sys.stdin,
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--mod_fp",
        help="File with config values to replace",
        type=argparse.FileType("r"),
        metavar="FILE",
    )
    group.add_argument(
        "--mod_str",
        help="Replacement YAML string (e.g. 'blast:{threads:4}')",
        type=str,
        metavar="STR",
    )

    args = parser.parse_args()
    config = ruamel.yaml.round_trip_load(args.config)

    if args.mod_fp:
        mods = ruamel.yaml.safe_load(args.mod_fp)
    else:
        mods = ruamel.yaml.safe_load(args.mod_str)
    print(mods)
    config = update_dict(config, mods)
    ruamel.yaml.round_trip_dump(config, sys.stdout)


if __name__ == "__main__":
    main()
