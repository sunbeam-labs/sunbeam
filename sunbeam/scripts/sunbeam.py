import argparse
import sys
from sunbeam import __version__
from sunbeam.scripts.config import main as Config
from sunbeam.scripts.extend import main as Extend
from sunbeam.scripts.init import main as Init
from sunbeam.scripts.run import main as Run


def main():
    parser = main_parser()
    args, remaining = parser.parse_known_args()

    if args.command == "run":
        Run(remaining)
    elif args.command == "init":
        Init(remaining)
    elif args.command == "config":
        Config(remaining)
    elif args.command == "extend":
        Extend(remaining)
    else:
        parser.print_help()
        sys.stderr.write("Unrecognized command.\n")


def main_parser():
    usage_str = "%(prog)s [-h/--help,-v/--version] <subcommand>"
    description_str = "Sunbeam subcommands: init, run, config, extend"

    parser = argparse.ArgumentParser(
        prog="sunbeam",
        usage=usage_str,
        description=description_str,
        epilog="For more help, see the docs at http://sunbeam.readthedocs.io.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
    )

    parser.add_argument("command", help=argparse.SUPPRESS, nargs="?")
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=__version__,
    )
    parser.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show this help message and exit.",
    )

    return parser
