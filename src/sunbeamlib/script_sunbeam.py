import argparse
import sys

import sunbeamlib
from sunbeamlib.script_run import main as Run
from sunbeamlib.script_init import main as Init
from sunbeamlib.script_config import main as Config
from sunbeamlib.script_list_samples import main as ListSamples
from sunbeamlib.script_extend import main as Extend


def main():
    usage_str = "%(prog)s [-h/--help,-v/--version] <subcommand>"
    description_str = (
        "subcommands:\n"
        "  init         \tCreate a new config file for a project using local or SRA data.\n"
        "  run          \tExecute the pipeline.\n"
        "  config       \tModify or update config files.\n"
        "  list_samples \tMake a list of samples from a directory.\n"
        "  extend       \tAdd an extension.\n"
    ).format(version=sunbeamlib.__version__)

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
        version=sunbeamlib.__version__,
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
