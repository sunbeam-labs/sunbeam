"""CLI for AI-driven workflow planning."""

import argparse
import sys

from sunbeam.ai import WorkflowPlanner


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Generate a bioinformatics workflow plan using an LLM.",
    )
    parser.add_argument("input_description", help="Description of the input data")
    parser.add_argument("output_description", help="Desired outputs")
    parser.add_argument(
        "--steps",
        default="",
        help="Comma separated list of allowed steps (optional)",
    )

    args = parser.parse_args(argv)

    planner = WorkflowPlanner()
    plan = planner.plan(args.input_description, args.output_description, args.steps)
    sys.stdout.write(plan + "\n")


if __name__ == "__main__":  # pragma: no cover
    main(sys.argv[1:])
