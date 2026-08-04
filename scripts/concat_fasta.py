#!/usr/bin/env python3
"""Concatenate FASTA inputs without changing headers or sequences."""

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("files", nargs="+")
    args = parser.parse_args()
    chunks = []
    for name in args.files:
        text = Path(name).read_text(encoding="utf-8-sig").strip()
        if text:
            chunks.append(text)
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("\n".join(chunks) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
