#!/usr/bin/env python3
"""Plot measured normalized CPM for pipeline-selected contact and loop candidates."""
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_selected(table, samples, destination, label):
    df = pd.read_csv(table, sep="\t")
    for sample in samples:
        if sample not in df.columns:
            raise ValueError(f"Missing sample-normalized column {sample} in {table}")
    df = df.sort_values("FDR", ascending=True).head(10)
    n = max(1, len(df))
    fig, axes = plt.subplots(n, 1, figsize=(10.5, min(28, 3.3*n+1.3)), squeeze=False)
    if df.empty:
        axes[0,0].text(0.5, 0.5, "No core candidates at configured thresholds", ha="center", va="center")
        axes[0,0].axis("off")
    for ax, row in zip(axes[:,0], df.itertuples(index=False)):
        record = row._asdict()
        values = [max(0.0, float(record[s])) for s in samples]
        ax.bar(range(len(samples)), values)
        ax.set_xticks(range(len(samples)), samples, rotation=25, ha="right")
        ax.set_ylabel("Normalized CPM")
        entity = record.get("contact_id", record.get("loop_id", "candidate"))
        desc = f"{entity}  {record.get('chrom1','')}:{record.get('start1','')}--{record.get('chrom2','')}:{record.get('start2','')}"
        ax.set_title(f"{desc}   FDR={float(record['FDR']):.3g}; logFC={float(record['logFC']):.2f}", fontsize=10)
    fig.suptitle(f"{label} - normalized counts by sample", fontsize=13)
    fig.text(0.03, 0.004, "Low-count focal candidates require manual checks; these are normalized CPM, not raw counts or direct proof of loop gains/losses.", fontsize=7.9)
    fig.tight_layout(rect=(0, .025, 1, .96))
    Path(destination).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--contacts", required=True)
    ap.add_argument("--loops", required=True)
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--out-contacts", required=True)
    ap.add_argument("--out-loops", required=True)
    args = ap.parse_args(argv)
    plot_selected(args.contacts, args.samples, args.out_contacts, "Differential contacts")
    plot_selected(args.loops, args.samples, args.out_loops, "Differential pooled-loop candidates")


if __name__ == "__main__":
    main()
