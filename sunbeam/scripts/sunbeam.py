import argparse
import sys
from sunbeam import __version__
from sunbeam.scripts.config import main as Config
from sunbeam.scripts.extend import main as Extend
from sunbeam.scripts.init import main as Init
from sunbeam.scripts.list_samples import main as ListSamples
from sunbeam.scripts.run import main as Run


def main():
    usage_str = "%(prog)s [-h/--help,-v/--version] <subcommand>"
    description_str = (
        "subcommands:\n"
        "  init         \tCreate a new config file for a project using local or SRA data.\n"
        "  run          \tExecute the pipeline.\n"
        "  config       \tModify or update config files.\n"
        "  list_samples \tMake a list of samples from a directory.\n"
        "  extend       \tAdd an extension.\n"
    ).format(version=__version__)

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

    args, remaining = parser.parse_known_args()

    if args.command == "run":
        Run(remaining)
    elif args.command == "init":
        Init(remaining)
    elif args.command == "config":
        Config(remaining)
    elif args.command == "list_samples":
        ListSamples(remaining)
    elif args.command == "extend":
        Extend(remaining)
    else:
        parser.print_help()
        sys.stderr.write("Unrecognized command.\n")
