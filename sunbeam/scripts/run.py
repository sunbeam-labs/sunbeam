import os
import sys
import argparse
import subprocess
from pathlib import Path
from sunbeam import __version__


def analyze_failure(log: str) -> None:
    """Use OpenAI to analyze failure logs."""
    try:
        from openai import OpenAI
    except ImportError:  # pragma: no cover - this is a soft dependency
        sys.stderr.write(
            "AI analysis requested, but the 'openai' package is not installed.\n"
        )
        return

    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key:
        sys.stderr.write("OPENAI_API_KEY not set; skipping AI analysis.\n")
        return

    try:
        client = OpenAI(api_key=api_key)
        resp = client.chat.completions.create(
            model="gpt-4.1-nano",
            messages=[
                {
                    "role": "system",
                    "content": "You diagnose errors from Sunbeam pipeline runs.",
                },
                {
                    "role": "user",
                    "content": f"Sunbeam run failed with the following output:\n{log}\n",
                },
            ],
            max_tokens=150,
        )
        sys.stderr.write("AI diagnosis:\n" + resp.choices[0].message.content + "\n")
    except Exception as exc:  # pragma: no cover - network errors are non-deterministic
        sys.stderr.write(f"AI analysis failed: {exc}\n")


def main(argv: list[str] = sys.argv):
    """CLI entry point for running Sunbeam."""
    parser = main_parser()
    args, remaining = parser.parse_known_args(argv)

    snakefile = Path(__file__).parent.parent / "workflow" / "Snakefile"
    if not snakefile.exists():
        sys.stderr.write(
            f"Error: could not find a Snakefile in directory '{snakefile}'\n"
        )
        sys.exit(1)

    conda_prefix = Path(__file__).parent.parent.parent / ".snakemake"

    conda_cmd = "conda" if not args.mamba else "mamba"

    if args.include and args.exclude:
        sys.stderr.write("Error: cannot use both --include and --exclude\n")
        sys.exit(1)

    os.environ["SUNBEAM_EXTS_INCLUDE"] = ""
    os.environ["SUNBEAM_EXTS_EXCLUDE"] = ""
    if args.include:
        os.environ["SUNBEAM_EXTS_INCLUDE"] = ", ".join(args.include)
    if args.exclude:
        os.environ["SUNBEAM_EXTS_EXCLUDE"] = ", ".join(args.exclude)

    if args.skip not in ["", "qc", "decontam"]:
        sys.stderr.write("Error: --skip must be either 'qc' or 'decontam'\n")
        sys.exit(1)

    os.environ["SUNBEAM_SKIP"] = args.skip

    os.environ["SUNBEAM_DOCKER_TAG"] = args.docker_tag

    # Extract the profile arg from the remaining args
    profile_parse = argparse.ArgumentParser(add_help=False)
    profile_parse.add_argument("--profile")
    profile_args, _ = profile_parse.parse_known_args(remaining)
    profile = profile_args.profile

    if not profile:
        sys.stderr.write(
            "Error: --profile is required. Please specify a profile to use.\n"
        )
        sys.exit(1)
    configfile = Path(profile) / "sunbeam_config.yml"

    snakemake_args = [
        "snakemake",
        "--snakefile",
        str(snakefile),
        "--conda-prefix",
        str(conda_prefix),
        "--conda-frontend",
        conda_cmd,
        "--configfile",
        str(configfile),
    ] + remaining
    sys.stderr.write("Running: " + " ".join(snakemake_args) + "\n")

    cmd = subprocess.run(snakemake_args, capture_output=True, text=True)
    if cmd.stdout:
        sys.stdout.write(cmd.stdout)
    if cmd.stderr:
        sys.stderr.write(cmd.stderr)
    if cmd.returncode != 0 and args.ai:
        analyze_failure(cmd.stdout + "\n" + cmd.stderr)

    sys.exit(cmd.returncode)


def main_parser():
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
        "-m",
        "--mamba",
        action="store_true",
        help="Use mamba instead of conda to create environments.",
    )
    parser.add_argument(
        "--include",
        action="append",
        default=[],
        help="Extension to include in run. Use multiple times for multiple extensions. Will exclude everything else.",
    )

    parser.add_argument(
        "--exclude",
        action="append",
        default=[],
        help="Extension to exclude from run. Use multiple times for multiple extensions. Use 'all' to exclude all extensions. Will include everything else.",
    )
    parser.add_argument(
        "--skip",
        default="",
        help="Workflow to skip. Either 'qc' to skip the quality control steps or 'decontam' to skip everything in sunbeam core (QC and decontamination).",
    )
    parser.add_argument(
        "--docker_tag",
        default=__version__,
        help="The tag to use when pulling docker images for the core pipeline environments, defaults to sunbeam's current version, a good alternative is 'latest' for the latest stable release",
    )
    parser.add_argument(
        "--ai",
        action="store_true",
        help="Use OpenAI to diagnose failures after the run",
    )

    return parser
