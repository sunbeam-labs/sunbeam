import argparse
from sunbeam.project import SunbeamConfig


def main(argv):
    parser = main_parser()
    args = parser.parse_args(argv)

    if not args.config_file:
        raise SystemExit("No config file specified. Use -h for help.")
    if args.update and args.modify:
        raise SystemExit(
            "Incompatible arguments: --update and --modify cannot be used together."
        )
    if args.update:
        sc = SunbeamConfig.from_file(args.config_file)
        sc.fill_missing()
        sc.to_file(args.config_file)
    if args.modify:
        sc = SunbeamConfig.from_file(args.config_file)
        sc.modify(args.modify)
        sc.to_file(args.config_file)
    if not args.update and not args.modify:
        raise SystemExit(
            "No action specified. Use --update to update the config file or "
            "--modify to modify it."
        )


def main_parser():
    parser = argparse.ArgumentParser(
        "sunbeam config",
        usage="%(prog)s [options]",
        description="Update or modify a Sunbeam config file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "config_file",
        type=str,
        help="Path to the Sunbeam config file to update or modify",
    )
    parser.add_argument(
        "-u",
        "--update",
        action="store_true",
        help="Update a config file to be compatible with the active version of Sunbeam and all installed extensions",
    )
    parser.add_argument(
        "-m",
        "--modify",
        type=str,
        help='Modify a config file with the specified changes (e.g. "qc: {minlen: 12}")',
    )

    return parser
