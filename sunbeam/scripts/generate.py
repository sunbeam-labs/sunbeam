import argparse
import re
import sys
from pathlib import Path

from sunbeam import EXTENSIONS_DIR, logger
from sunbeam.ai.rule_creator import create_rules_from_prompt


def main(argv=sys.argv):
    parser = main_parser()
    args = parser.parse_args(argv)

    ruleset_name = args.name
    prompt = args.prompt

    ext_dir = EXTENSIONS_DIR() / f"sbx_{ruleset_name}"
    try:
        ext_dir.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        raise SystemExit(f"Extension directory {ext_dir} already exists")

    rules_path = ext_dir / f"sbx_{ruleset_name}.smk"

    # Default context files to include
    context_files = [
        Path(__file__).parent.parent / "workflow" / "Snakefile",
        Path(__file__).parent.parent / "workflow" / "rules" / "qc.smk",
        Path(__file__).parent.parent / "workflow" / "rules" / "decontam.smk",
    ]

    result = create_rules_from_prompt(
        ruleset_name, prompt, context_files=context_files, write_to=rules_path
    )

    logger.info(f"Created extension scaffold at {ext_dir}")


def main_parser():
    parser = argparse.ArgumentParser(
        "sunbeam generate",
        description="Generate a Sunbeam extension from a natural language prompt.",
    )
    parser.add_argument(
        "--name", "-n", required=True, help="Name of the ruleset to create"
    )
    parser.add_argument(
        "--prompt", "-p", required=True, help="Description of rules to generate"
    )
    return parser
