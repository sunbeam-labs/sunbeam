import sys
import argparse
import subprocess
from sunbeam import EXTENSIONS_DIR, logger


def main(argv=sys.argv):
    parser = main_parser()
    args = parser.parse_args(argv)

    if args.github_url.startswith("http"):
        gh_url = args.github_url
        ext_name = gh_url.split("/")[-1].replace(".git", "")
        logger.info(
            f"GitHub URL provided. Cloning from the specified URL: {gh_url}. Parsed extension name: {ext_name}"
        )
    else:
        ext_name = args.github_url
        gh_url = "https://github.com/sunbeam-labs/" + ext_name + ".git"
        logger.info(
            f"Extension name provided. Constructed GitHub URL: {gh_url}. Using extension name: {ext_name}"
        )

    git_clone_args = ["git", "clone", gh_url, str(EXTENSIONS_DIR() / ext_name)]
    cmd = subprocess.run(git_clone_args)


def main_parser():
    parser = argparse.ArgumentParser(
        "sunbeam extend",
        usage="%(prog)s github_url",
        description="Installs a Sunbeam extension from the given name or GitHub URL.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "github_url",
        type=str,
        help="Name of extension or GitHub URL.",
    )

    return parser
