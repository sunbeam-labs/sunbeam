import os
import sys
import argparse
import subprocess
from pathlib import Path


def main(argv=sys.argv):

    epilog_str = (
        "You can pass further arguments to Snakemake after --, e.g:\n"
        "    $ sunbeam run -- --cores 12\n"
        "See http://snakemake.readthedocs.io for more information.\n"
        " "
    )

    parser = argparse.ArgumentParser(
        "sunbeam run",
        usage="%(prog)s [options] -- <snakemake options>",
        description="Executes the Sunbeam pipeline by calling Snakemake.",
        epilog=epilog_str,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-s",
        "--sunbeam_dir",
        default=os.getenv("SUNBEAM_DIR", os.getcwd()),
        help="Path to Sunbeam installation",
    )
    parser.add_argument(
        "--target_list", nargs="+", default=[], help="List of sunbeam targets"
    )

    # The remaining args (after --) are passed to Snakemake
    args, remaining = parser.parse_known_args(argv)

    snakefile = Path(args.sunbeam_dir) / "Snakefile"
    if not snakefile.exists():
        sys.stderr.write(
            f"Error: could not find a Snakefile in directory '{args.sunbeam_dir}'\n"
        )
        sys.exit(1)

    conda_prefix = Path(args.sunbeam_dir) / ".snakemake"

    cmds = list()
    if args.target_list == []:
        args.target_list = [""]

    for target in args.target_list:
        if target:
            print(f"Running sunbeam on target: {target}")

        # Including target when it's en empty string breaks stuff so the extra
        # list comp avoids that
        snakemake_args = [
            arg
            for arg in [
                "snakemake",
                "--snakefile",
                str(snakefile),
                "--conda-prefix",
                str(conda_prefix),
                target,
            ]
            if arg
        ] + remaining
        print("Running: " + " ".join(snakemake_args))

        cmds.append(subprocess.run(snakemake_args))

    sys.exit(cmds[0].returncode)
