import argparse
import contextlib
import datetime
import logging
import os
import sys
from pathlib import Path
from snakemake.cli import main as snakemake_main
from sunbeam import __version__
from sunbeam.logging import get_pipeline_logger, StreamToLogger


def analyze_run(log: str, logger: logging.Logger, ai: bool) -> None:
    """Analyze the run log and provide insights or suggestions."""
    # We could do some rule-based analysis here but I'd rather lean into the AI features and see how far they can take us
    if ai:
        try:
            from openai import OpenAI
        except ImportError:  # pragma: no cover - this is a soft dependency
            logger.error(
                "AI analysis requested, but the 'openai' package is not installed.\n"
            )
            return

        api_key = os.getenv("OPENAI_API_KEY")
        if not api_key:
            logger.error("OPENAI_API_KEY not set; skipping AI analysis.\n")
            return

        try:
            client = OpenAI(api_key=api_key)
            resp = client.chat.completions.create(
                model="gpt-4.1-nano",
                messages=[
                    {
                        "role": "system",
                        "content": "You diagnose errors from Sunbeam pipeline runs. If there are problems, suggest possible causes and solutions. Keep the answer short and sweet. If there are relevant file paths for debugging (like log files), mention them.",
                    },
                    {
                        "role": "user",
                        "content": f"Sunbeam ran with the following output:\n{log}\n",
                    },
                ],
                max_tokens=1500,
            )
            logger.info(
                "\n\nAI diagnosis:\n"
                + resp.choices[0].message.content
                + "\nCheck out the Sunbeam documentation (https://sunbeam.readthedocs.io/en/stable/) and the GitHub issues page (https://github.com/sunbeam-labs/sunbeam/issues) for more information or to open a new issue.\n"
            )
        except (
            Exception
        ) as exc:  # pragma: no cover - network errors are non-deterministic
            logger.error(f"AI analysis failed: {exc}\n")


def main(argv: list[str] = sys.argv):
    """CLI entry point for running Sunbeam."""
    parser = main_parser()
    args, remaining = parser.parse_known_args(argv)

    profile = Path(args.profile).resolve()

    log_file = args.log_file
    if not log_file:
        log_file = (
            profile
            / f"sunbeam_{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log"
        )

    # From here on everything is considered part of the "pipeline"
    # This means all logs are handled by the pipeline logger (or pipeline extension loggers)
    # You could argue it would make more sense to start this at the actual snakemake call
    # but this way we can log some relevant setup information that might be useful on post-mortem analysis
    logger = get_pipeline_logger(log_file)
    logger.debug("Sunbeam pipeline logger initialized.")

    snakefile = Path(__file__).parent.parent / "workflow" / "Snakefile"
    if not snakefile.exists():
        logger.error(f"Could not find a Snakefile in directory '{snakefile}'\n")
        sys.exit(1)

    conda_prefix = Path(__file__).parent.parent.parent / ".snakemake"

    conda_cmd = "conda" if not args.mamba else "mamba"

    if args.include and args.exclude:
        logger.error("Cannot use both --include and --exclude\n")
        sys.exit(1)

    os.environ["SUNBEAM_EXTS_INCLUDE"] = ""
    os.environ["SUNBEAM_EXTS_EXCLUDE"] = ""
    if args.include:
        os.environ["SUNBEAM_EXTS_INCLUDE"] = ", ".join(args.include)
    if args.exclude:
        os.environ["SUNBEAM_EXTS_EXCLUDE"] = ", ".join(args.exclude)

    if args.skip not in ["", "qc", "decontam"]:
        logger.error("The value of --skip must be either 'qc' or 'decontam'\n")
        sys.exit(1)

    os.environ["SUNBEAM_SKIP"] = args.skip

    os.environ["SUNBEAM_DOCKER_TAG"] = args.docker_tag

    configfile = Path(profile) / "sunbeam_config.yml"

    snakemake_args = [
        "--snakefile",
        str(snakefile),
        "--profile",
        str(profile),
        "--configfile",
        str(configfile),
        "--conda-prefix",
        str(conda_prefix),
        "--conda-frontend",
        conda_cmd,
    ] + remaining
    logger.info("Running: " + " ".join(snakemake_args))

    try:
        stream_logger = StreamToLogger(logger, level=logging.INFO)

        with contextlib.redirect_stderr(stream_logger):
            snakemake_main(snakemake_args)
    finally:
        print(log_file)
        print(log_file.exists())
        with open(log_file, "r") as f:
            analyze_run(f.read(), logger, args.ai)


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
        "--profile",
        required=True,
        help="The Snakemake profile to use for running the pipeline. This should be a directory containing a 'config.yaml' file.",
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
    parser.add_argument(
        "--log_file",
        default=None,
        help="Path to a file where the pipeline log will be written. If not specified, logs will be written to `sunbeam_<current_datetime>.log` under the project directory.",
    )

    return parser
