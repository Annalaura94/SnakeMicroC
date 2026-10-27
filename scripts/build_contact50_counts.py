#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import cooler
import numpy as np
import pandas as pd


def parse_csv(x):
    return [v for v in x.split(",") if v]


def ordered_candidate_keys(sample_dicts):
    """Return genomic pixel coordinates in a stable, reproducible order.

    Iterating over a set would assign different contact IDs across runs,
    even when the underlying contact matrix and counts were identical.
    """
    return sorted(set().union(*(d.keys() for d in sample_dicts.values())))


def chromosome_names(view_path):
    x = pd.read_csv(view_path, sep="\t", header=None)
    return x.iloc[:, 0].astype(str).tolist()


def sparse_upper_counts(clr, chrom, res, min_dist, max_dist):
    m = clr.matrix(
        balance=False,
        sparse=True
    ).fetch(chrom).tocoo()

    keep = m.row < m.col

    row = m.row[keep].astype(np.int64)
    col = m.col[keep].astype(np.int64)
    data = m.data[keep].astype(np.int64)

    dist = (col - row) * res

    keep = (
        (dist >= min_dist) &
        (dist <= max_dist) &
        (data > 0)
    )

    return {
        (int(i), int(j)): int(v)
        for i, j, v in zip(
            row[keep],
            col[keep],
            data[keep]
        )
    }


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--samples", required=True)
    ap.add_argument("--resolution", type=int, required=True)
    ap.add_argument("--min-distance-bp", type=int, required=True)
    ap.add_argument("--max-distance-bp", type=int, required=True)
    ap.add_argument("--prefilter-min-total-count", type=int, required=True)
    ap.add_argument("--prefilter-min-nonzero-samples", type=int, required=True)
    ap.add_argument("--mcool-pattern", required=True)
    ap.add_argument("--autosomes-view", required=True)
    ap.add_argument("--outdir", required=True)

    args = ap.parse_args()

    samples = parse_csv(args.samples)

    if len(samples) < 2:
        raise ValueError("At least two samples are required")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    chroms = chromosome_names(
        Path(args.autosomes_view)
    )

    clrs = {
        s: cooler.Cooler(
            f"{args.mcool_pattern.format(sample=s)}"
            f"::/resolutions/{args.resolution}"
        )
        for s in samples
    }

    out_path = (
        outdir /
        "contact50_condition_blind_counts.tsv.gz"
    )

    summary_rows = []
    first = True
    contact_counter = 0

    for chrom in chroms:
        bins = (
            clrs[samples[0]]
            .bins()
            .fetch(chrom)
            .reset_index(drop=True)
        )

        sample_dicts = {
            s: sparse_upper_counts(
                clrs[s],
                chrom,
                args.resolution,
                args.min_distance_bp,
                args.max_distance_bp
            )
            for s in samples
        }

        keys = ordered_candidate_keys(sample_dicts)

        n_union = len(keys)
        rows = []

        for i, j in keys:
            vals = [
                sample_dicts[s].get((i, j), 0)
                for s in samples
            ]

            total = int(sum(vals))
            nnz = int(sum(v > 0 for v in vals))

            if (
                total <
                args.prefilter_min_total_count
                or
                nnz <
                args.prefilter_min_nonzero_samples
            ):
                continue

            contact_counter += 1

            rec = {
                "contact_id":
                    f"C50_{contact_counter:08d}",
                "chrom1": chrom,
                "start1": int(
                    bins.iloc[i]["start"]
                ),
                "end1": int(
                    bins.iloc[i]["end"]
                ),
                "chrom2": chrom,
                "start2": int(
                    bins.iloc[j]["start"]
                ),
                "end2": int(
                    bins.iloc[j]["end"]
                ),
                "distance_bp":
                    int((j - i) * args.resolution),
                "total_raw": total,
                "nonzero_samples": nnz,
            }

            for s, v in zip(samples, vals):
                rec[s] = int(v)

            rows.append(rec)

        d = pd.DataFrame(rows)

        if len(d):
            d.to_csv(
                out_path,
                sep="\t",
                index=False,
                compression="gzip",
                mode="wt" if first else "at",
                header=first,
            )

            first = False

        summary_rows.append({
            "chrom": chrom,
            "union_nonzero_pixels": n_union,
            "retained_condition_blind_prefilter":
                len(d),
        })

        print(
            f"{chrom}: union={n_union:,}; "
            f"retained={len(d):,}"
        )

    if first:
        raise RuntimeError(
            "No contacts passed the "
            "condition-blind prefilter"
        )

    summary = pd.DataFrame(summary_rows)

    summary.loc[len(summary)] = {
        "chrom": "TOTAL",
        "union_nonzero_pixels":
            int(
                summary[
                    "union_nonzero_pixels"
                ].sum()
            ),
        "retained_condition_blind_prefilter":
            int(
                summary[
                    "retained_condition_blind_prefilter"
                ].sum()
            ),
    }

    summary.to_csv(
        outdir /
        "contact50_prefilter_summary.tsv",
        sep="\t",
        index=False
    )

    print(f"Written {out_path}")


if __name__ == "__main__":
    main()
