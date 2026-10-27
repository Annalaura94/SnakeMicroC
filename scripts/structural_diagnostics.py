#!/usr/bin/env python3
"""Replicate-concordant descriptive E1 switches and insulation boundary calls.

These screening calls are NOT differential significance tests. The E1 sign is
assumed already oriented by cooltools using GC. Only bin positions passing
threshold in every listed sample of each group are called concordant switches.
"""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def strict_switches(df, ctrl, test, threshold):
    needed = {"chrom", "start", "end", *ctrl, *test}
    if not needed.issubset(df.columns):
        raise ValueError(f"E1 summary missing columns: {sorted(needed - set(df.columns))}")
    c = df[ctrl].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    t = df[test].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    a_b = np.isfinite(c).all(axis=1) & np.isfinite(t).all(axis=1)
    a_to_b = a_b & (c >= threshold).all(axis=1) & (t <= -threshold).all(axis=1)
    b_to_a = a_b & (c <= -threshold).all(axis=1) & (t >= threshold).all(axis=1)
    out = df.loc[a_to_b | b_to_a, ["chrom", "start", "end", *ctrl, *test]].copy()
    out["switch"] = np.where(a_to_b[a_to_b | b_to_a], "A_to_B", "B_to_A")
    out["mean_control_e1"] = np.mean(c[a_to_b | b_to_a], axis=1)
    out["mean_test_e1"] = np.mean(t[a_to_b | b_to_a], axis=1)
    out["delta_e1"] = out["mean_test_e1"] - out["mean_control_e1"]
    return out.sort_values("delta_e1", key=abs, ascending=False)


def read_boundary_points(path, sample, window):
    df = pd.read_csv(path, sep="\t")
    name = f"is_boundary_{window}"
    if name not in df.columns:
        candidates = [c for c in df if c.startswith("is_boundary_") and str(window) in c]
        if len(candidates) != 1:
            raise ValueError(f"{path}: no unambiguous boundary column for {window}: {candidates}")
        name = candidates[0]
    observed = df[name].astype(str).str.strip().str.lower().isin(("true", "1"))
    rows = df.loc[observed, ["chrom", "start", "end"]].copy()
    rows["sample"] = sample
    rows["center"] = ((rows["start"] + rows["end"]) // 2).astype(int)
    return rows


def cluster_boundaries(boundaries, samples, ctrl, test, tolerance):
    columns = ["chrom", "start", "end", "support_control", "support_test",
               "support_total", "event", "sample_support"]
    if not boundaries:
        return pd.DataFrame(columns=columns)
    points = pd.concat(boundaries, ignore_index=True)
    if points.empty:
        return pd.DataFrame(columns=columns)
    result = []
    for chrom, tab in points.groupby("chrom", sort=True):
        tab = tab.sort_values("center")
        clusters = []
        for row in tab.itertuples(index=False):
            pos = int(row.center)
            if not clusters or pos - clusters[-1][0] > tolerance:
                clusters.append((pos, [(pos, row.sample)]))
            else:
                clusters[-1][1].append((pos, row.sample))
        for _, members in clusters:
            center = int(np.median([x[0] for x in members]))
            present = {x[1] for x in members}
            ncontrol = len(set(ctrl) & present)
            ntest = len(set(test) & present)
            if ncontrol == 0 and ntest == len(test):
                event = "gain_in_test"
            elif ncontrol == len(ctrl) and ntest == 0:
                event = "loss_in_test"
            elif ncontrol == len(ctrl) and ntest == len(test):
                event = "shared"
            else:
                event = "mixed_or_single_sample"
            result.append(dict(chrom=chrom, start=center, end=center+1,
                               support_control=ncontrol, support_test=ntest,
                               support_total=len(present), event=event,
                               sample_support=",".join(sorted(present))))
    return pd.DataFrame(result, columns=columns)


def plot_counts(counts, labels, title, path, note):
    fig, ax = plt.subplots(figsize=(8.5, 5.1))
    vals = [int(counts.get(x, 0)) for x in labels]
    ax.bar(range(len(vals)), vals)
    ax.set_xticks(range(len(labels)), [x.replace("_", " ") for x in labels], rotation=20, ha="right")
    ax.set_ylabel("Bins / boundary clusters")
    ax.set_title(title)
    for i, v in enumerate(vals):
        ax.text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=9)
    fig.text(0.03, 0.015, note, fontsize=8, color="#444")
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(path)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--e1-summary", required=True)
    ap.add_argument("--insulation", nargs="+", required=True)
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--control-samples", nargs="+", required=True)
    ap.add_argument("--test-samples", nargs="+", required=True)
    ap.add_argument("--e1-threshold", type=float, required=True)
    ap.add_argument("--insulation-window", type=int, required=True)
    ap.add_argument("--tolerance", type=int, required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args(argv)
    if len(args.samples) != len(args.insulation):
        ap.error("One insulation track per sample required")
    if len(args.control_samples) < 2 or len(args.test_samples) < 2:
        ap.error("At least two samples per comparison group required")
    if set(args.control_samples + args.test_samples) - set(args.samples):
        ap.error("Control/test group contains unknown samples")
    if args.e1_threshold <= 0 or args.tolerance < 0:
        ap.error("E1 threshold must be >0 and tolerance >=0")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    e1 = pd.read_csv(args.e1_summary, sep="\t")
    switches = strict_switches(e1, args.control_samples, args.test_samples, args.e1_threshold)
    switches.to_csv(outdir / "compartment_switches.tsv.gz", sep="\t", index=False, compression="gzip")
    sw_counts = {k: int((switches["switch"] == k).sum()) for k in ("A_to_B", "B_to_A")}
    sw_summary = pd.DataFrame([{"metric": "evaluated_bins", "value": len(e1)},
                               {"metric": "strict_A_to_B", "value": sw_counts["A_to_B"]},
                               {"metric": "strict_B_to_A", "value": sw_counts["B_to_A"]},
                               {"metric": "absolute_e1_threshold", "value": args.e1_threshold}])
    sw_summary.to_csv(outdir / "compartment_switch_summary.tsv", sep="\t", index=False)
    plot_counts(sw_counts, ["A_to_B", "B_to_A"], "Strict replicate-concordant E1 sign switches",
                outdir / "compartment_switches.pdf", "Descriptive calls in all samples; not statistical discoveries.")
    points = [read_boundary_points(path, s, args.insulation_window)
              for path, s in zip(args.insulation, args.samples)]
    events = cluster_boundaries(points, args.samples, args.control_samples, args.test_samples, args.tolerance)
    events.to_csv(outdir / "boundary_events.tsv.gz", sep="\t", index=False, compression="gzip")
    keys = ("gain_in_test", "loss_in_test", "shared", "mixed_or_single_sample")
    counts = {k: int((events["event"] == k).sum()) for k in keys}
    bd_summary = pd.DataFrame([{"metric": "boundary_clusters", "value": len(events)},
                               *({"metric": k, "value": v} for k, v in counts.items()),
                               {"metric": "cluster_tolerance_bp", "value": args.tolerance}])
    bd_summary.to_csv(outdir / "boundary_events_summary.tsv", sep="\t", index=False)
    plot_counts(counts, keys, "Replicate-concordant insulation boundary events",
                outdir / "boundary_events.pdf", "Boundary presence is descriptive; no differential FDR is calculated.")
    print(f"E1 switches {sw_counts}; boundary events {counts}")


if __name__ == "__main__":
    main()
