import os
import sys
import argparse
import subprocess
from pathlib import Path

from sunbeamlib import __version__


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
        "-m",
        "--mamba",
        action="store_true",
        help="Use mamba instead of conda to create environments",
    )
    parser.add_argument(
        "--target_list",
        nargs="+",
        default=[],
        help="List of sunbeam targets (DEPRECATED)",
    )
    parser.add_argument(
        "--include",
        nargs="+",
        default=[],
        help="List of extensions to include in run",
    )
    parser.add_argument(
        "--exclude",
        nargs="+",
        default=[],
        help="List of extensions to exclude from run, use 'all' to exclude all extensions",
    )
    parser.add_argument(
        "--docker_tag",
        default=__version__,
        help="The tag to use when pulling docker images for the core pipeline environments, defaults to sunbeam's current version ($SUNBEAM_VER), a good alternative is 'latest' for the latest stable release",
    )

    # The remaining args (after --) are passed to Snakemake
    args, remaining = parser.parse_known_args(argv)

    snakefile = Path(args.sunbeam_dir) / "workflow" / "Snakefile"
    if not snakefile.exists():
        sys.stderr.write(
            f"Error: could not find a Snakefile in directory '{str(Path(args.sunbeam_dir) / 'workflow')}'\n"
        )
        sys.exit(1)

    conda_prefix = Path(args.sunbeam_dir) / ".snakemake"

    conda_cmd = "conda" if not args.mamba else "mamba"

    if args.target_list:
        sys.stderr.write(
            "Warning: passing targets to '--target_list' is deprecated. "
            "Please use 'sunbeam run <opts> target1 target2 target3' instead.\n"
        )

    if args.include and args.exclude:
        sys.stderr.write("Error: cannot use both --include and --exclude\n")
        sys.exit(1)

    os.environ["SUNBEAM_EXTS_INCLUDE"] = ""
    os.environ["SUNBEAM_EXTS_EXCLUDE"] = ""
    if args.include:
        os.environ["SUNBEAM_EXTS_INCLUDE"] = ", ".join(args.include)
    if args.exclude:
        os.environ["SUNBEAM_EXTS_EXCLUDE"] = ", ".join(args.exclude)

    os.environ["SUNBEAM_DOCKER_TAG"] = args.docker_tag

    snakemake_args = (
        [
            "snakemake",
            "--snakefile",
            str(snakefile),
            "--conda-prefix",
            str(conda_prefix),
            "--conda-frontend",
            conda_cmd,
        ]
        + remaining
        + args.target_list
    )
    sys.stderr.write("Running: " + " ".join(snakemake_args))

    cmd = subprocess.run(snakemake_args)

    sys.exit(cmd.returncode)
