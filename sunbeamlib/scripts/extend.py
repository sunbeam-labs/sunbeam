import os
import sys
import argparse
import subprocess
from pathlib import Path

def main(argv=sys.argv):

    parser = argparse.ArgumentParser(
        "sunbeam extend",
        usage="%(prog)s github_url",
        description="Installs a Sunbeam extension from the given GitHub URL.",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        "github_url", type = str,
        help="URL to Sunbeam extension on GitHub (e.g. https://github.com/sunbeam-labs/sbx_report")

    parser.add_argument(
        "-s", "--sunbeam_dir", default=os.getenv("SUNBEAM_DIR", os.getcwd()),
        help="Path to Sunbeam installation")

    args = parser.parse_args(argv)
    
    extensions_dir = Path(args.sunbeam_dir)/"extensions"
    if not extensions_dir.exists():
        sys.stderr.write(
            "Error: could not find an extensions directory in '{}'\n".format(
                args.sunbeam_dir))
        sys.exit(1)

    extension_name = args.github_url.split("/")[-1]
    if extension_name.endswith(".git"):
        extension_name = extension_name[:-4]

    git_clone_args = ["git","clone",args.github_url,str(extensions_dir/extension_name)]

    print("Installing: "+extension_name+" from "+args.github_url))

    cmd = subprocess.run(git_clone_args)
    
    sys.exit(cmd.returncode)

