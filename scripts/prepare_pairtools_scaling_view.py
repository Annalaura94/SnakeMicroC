#!/usr/bin/env python3
"""Write a named-column view for pairtools scaling without changing BED4 views.

The canonical main view (for cooltools and cooler) is headerless BED4:
    chr1\t0\t248956422\tchr1

Pairtools 1.1.3 scaling reads a headered TSV with columns:
    chrom\tstart\tend\tname

Do not modify the canonical view because other tools rely on headerless BED4.
"""
from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path

HEADER = ("chrom", "start", "end", "name")


def make_view(input_path: str | Path, output_path: str | Path) -> int:
    source = Path(input_path)
    dest = Path(output_path)
    if source.resolve() == dest.resolve():
        raise ValueError("Input and output view paths must be different")

    records = []
    names = set()
    with source.open(newline="") as handle:
        for line_no, fields in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not fields or (len(fields) == 1 and not fields[0].strip()):
                continue
            if fields[0].startswith("#"):
                continue
            if tuple(fields) == HEADER and not records:
                continue
            if len(fields) != 4:
                raise ValueError(f"{source}:{line_no}: expected BED4, got {len(fields)} fields")
            chrom, start_s, end_s, name = [f.strip() for f in fields]
            try:
                start, end = int(start_s), int(end_s)
            except ValueError as exc:
                raise ValueError(f"{source}:{line_no}: start and end must be integers") from exc
            if not chrom or not name or start < 0 or end <= start:
                raise ValueError(f"{source}:{line_no}: invalid chromosome, name or coordinates")
            if name in names:
                raise ValueError(f"{source}:{line_no}: duplicated view region name {name!r}")
            names.add(name)
            records.append((chrom, start, end, name))

    if not records:
        raise ValueError(f"{source}: no valid regions")

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp_name = None
    try:
        with tempfile.NamedTemporaryFile(
            "w", newline="", dir=dest.parent, prefix=f".{dest.name}.", delete=False
        ) as handle:
            tmp_name = handle.name
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(HEADER)
            writer.writerows(records)
        os.replace(tmp_name, dest)
    finally:
        if tmp_name and os.path.exists(tmp_name):
            os.unlink(tmp_name)
    return len(records)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    n = make_view(args.input, args.output)
    print(f"Wrote {n} named view regions to {args.output}")


if __name__ == "__main__":
    main()
