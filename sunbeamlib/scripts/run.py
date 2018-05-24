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
        " ")

    parser = argparse.ArgumentParser(
        "sunbeam run",
        usage="%(prog)s [options] -- <snakemake options>",
        description="Executes the Sunbeam pipeline by calling Snakemake.",
        epilog=epilog_str,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        "-s", "--sunbeam_dir", default=os.getenv("SUNBEAM_DIR", os.getcwd()),
        help="Path to Sunbeam installation")

    # The remaining args (after --) are passed to Snakemake
    args, remaining = parser.parse_known_args(argv)

    snakefile = Path(args.sunbeam_dir)/"Snakefile"
    if not snakefile.exists():
        sys.stderr.write(
            "Error: could not find a Snakefile in directory '{}'\n".format(
                args.sunbeam_dir))
        sys.exit(1)
        
    snakemake_args = ['snakemake', '--snakefile', str(snakefile)] + remaining
    print("Running: "+" ".join(snakemake_args))

    cmd = subprocess.run(snakemake_args)
    
    sys.exit(cmd.returncode)
    
