import os
import sys
import argparse
from pathlib import Path

def main():

    parser = argparse.ArgumentParser(
        description="Finds target files; breaks if not found.")
    parser.add_argument(
        "--prefix", type=Path, help="Prefix to add to paths in target",
        default="")
    parser.add_argument(
        "targets", type=argparse.FileType('r'))

    args = parser.parse_args()
    for line in args.targets:
        if not line.strip():
            continue
        target = args.prefix/Path(line.strip())
        if not target.exists():
            raise SystemExit("Target '{}' not found".format(target))
        else:
            print("Found target '{}'".format(target))


if __name__ == "__main__":
    main()
