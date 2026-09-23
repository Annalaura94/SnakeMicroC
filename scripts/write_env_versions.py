#!/usr/bin/env python3

import argparse
import importlib.metadata
import platform
import subprocess
import pandas as pd


def command(cmd):
    try:
        return subprocess.check_output(
            cmd,
            shell=True,
            stderr=subprocess.STDOUT,
            text=True,
        ).strip()
    except Exception as exc:
        return f"ERROR: {exc}"


def pkg(name):
    try:
        return importlib.metadata.version(name)
    except Exception:
        return "not_available"


ap = argparse.ArgumentParser()
ap.add_argument("--mode", choices=["processing", "downstream"], required=True)
ap.add_argument("--output", required=True)
args = ap.parse_args()

rows = []
rows.append((args.mode, "python", platform.python_version()))

if args.mode == "processing":
    rows.extend([
        ("processing", "bwa", command("bwa 2>&1 | grep -m1 Version")),
        ("processing", "samtools", command("samtools --version | head -1")),
        ("processing", "pairtools", pkg("pairtools")),
        ("processing", "pysam", pkg("pysam")),
    ])

elif args.mode == "downstream":
    rows.extend([
        ("downstream", "cooler", pkg("cooler")),
        ("downstream", "cooltools", pkg("cooltools")),
        ("downstream", "bioframe", pkg("bioframe")),
        ("downstream", "scikit-image", pkg("scikit-image")),
        ("downstream", "hicrep", pkg("hicrep")),
        ("downstream", "numpy", pkg("numpy")),
        ("downstream", "pandas", pkg("pandas")),
    ])

pd.DataFrame(
    rows,
    columns=["environment", "software", "version"],
).to_csv(
    args.output,
    sep="\t",
    index=False,
)
