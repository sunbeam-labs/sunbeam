import argparse
import subprocess
import sunbeamlib
from sunbeamlib.scripts.run import main as Run

def main():

    usage_str = "%(prog)s [-h/--help,-v/--version] <subcommand>"
    description_str = (
        "subcommands:\n"
        "  init   \tCreate a new config file for a project.\n"
        "  run    \tExecute the pipeline.\n"
        "  config \tModify or update config files.\n"
    ).format(version=sunbeamlib.__version__)

    parser = argparse.ArgumentParser(
        prog="sunbeam",
        usage=usage_str,
        description=description_str,
        epilog="For more help, see the docs at http://sunbeam.readthedocs.io.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False
    )

    parser.add_argument("command", help=argparse.SUPPRESS, nargs="?")
    parser.add_argument(
        "-v", "--version", action="version",
        version="%(prog)s v{}".format(sunbeamlib.__version__))

    args, remaining = parser.parse_known_args()

    if args.command == "run":
        Run(remaining)
    elif args.command is None:
        parser.print_help()
        
