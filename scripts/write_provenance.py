#!/usr/bin/env python3

from pathlib import Path
import json
import pandas as pd
import yaml


def flatten(obj, prefix=""):
    rows = []
    if isinstance(obj, dict):
        for key, value in obj.items():
            new_key = f"{prefix}.{key}" if prefix else str(key)
            rows.extend(flatten(value, new_key))
    elif isinstance(obj, list):
        rows.append((prefix, json.dumps(obj)))
    else:
        rows.append((prefix, obj))
    return rows


config_out = Path(snakemake.output.config)
versions_out = Path(snakemake.output.versions)
summary_out = Path(snakemake.output.summary)

config_out.parent.mkdir(parents=True, exist_ok=True)

with config_out.open("w") as handle:
    yaml.safe_dump(dict(snakemake.config), handle, sort_keys=False)

version_tables = [
    pd.read_csv(snakemake.input.processing, sep="\t"),
    pd.read_csv(snakemake.input.downstream, sep="\t"),
    pd.read_csv(snakemake.input.r, sep="\t"),
]

pd.concat(version_tables, ignore_index=True).to_csv(
    versions_out,
    sep="\t",
    index=False,
)

pd.DataFrame(
    flatten(dict(snakemake.config)),
    columns=["parameter", "value"],
).to_csv(
    summary_out,
    sep="\t",
    index=False,
)
