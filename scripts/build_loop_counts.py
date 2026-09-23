#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import cooler
import numpy as np
import pandas as pd


def parse_csv(x):
    return [v for v in x.split(",") if v]


def window_sum(mat, i, j, r):
    i0 = max(0, i-r)
    i1 = min(mat.shape[0], i+r+1)

    j0 = max(0, j-r)
    j1 = min(mat.shape[1], j+r+1)

    return float(
        mat[i0:i1, j0:j1].sum()
    )


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--samples", required=True)
    ap.add_argument("--resolution", type=int, required=True)
    ap.add_argument("--min-distance-bp", type=int, required=True)
    ap.add_argument("--max-distance-bp", type=int, required=True)
    ap.add_argument("--quant-radius-bins", type=int, required=True)

    ap.add_argument("--dots", required=True)
    ap.add_argument("--cool-pattern", required=True)
    ap.add_argument("--outdir", required=True)

    args = ap.parse_args()

    samples = parse_csv(args.samples)

    if len(samples) < 2:
        raise ValueError(
            "At least two samples are required"
        )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    dots = pd.read_csv(
        args.dots,
        sep="\t"
    )

    req = [
        "chrom1",
        "start1",
        "end1",
        "chrom2",
        "start2",
        "end2",
    ]

    missing = [
        x for x in req
        if x not in dots.columns
    ]

    if missing:
        raise ValueError(
            f"Dot table missing {missing}; "
            f"columns={list(dots.columns)}"
        )

    dots = dots[
        dots.chrom1.eq(dots.chrom2)
    ].copy()

    dots["distance_bp"] = (
        dots.start2 - dots.start1
    ).abs()

    dots = dots[
        dots.distance_bp.between(
            args.min_distance_bp,
            args.max_distance_bp
        )
    ].copy()

    dots = (
        dots
        .sort_values(
            ["chrom1", "start1", "start2"]
        )
        .reset_index(drop=True)
    )

    dots["loop_id"] = [
        f"LOOP_{i:06d}"
        for i in range(
            1,
            len(dots) + 1
        )
    ]

    cools = {
        s: cooler.Cooler(
            args.cool_pattern.format(
                sample=s,
                resolution=args.resolution
            )
        )
        for s in samples
    }

    outrows = []

    for chrom, d in dots.groupby(
        "chrom1",
        sort=False
    ):
        cache = {}

        for s, clr in cools.items():
            bins = (
                clr.bins()
                .fetch(chrom)
                .reset_index(drop=True)
            )

            start_to_idx = {
                int(x): i
                for i, x in enumerate(
                    bins.start
                )
            }

            raw = (
                clr.matrix(
                    balance=False,
                    sparse=True
                )
                .fetch(chrom)
                .tocsr()
            )

            bal = (
                clr.matrix(
                    balance=True,
                    sparse=True
                )
                .fetch(chrom)
                .tocsr()
            )

            cache[s] = (
                start_to_idx,
                raw,
                bal
            )

        for rr in d.itertuples(
            index=False
        ):
            rec = {
                k: getattr(rr, k)
                for k in req
            }

            rec["loop_id"] = rr.loop_id
            rec["distance_bp"] = int(
                rr.distance_bp
            )

            for s in samples:
                idx, raw, bal = cache[s]

                if (
                    int(rr.start1) not in idx
                    or
                    int(rr.start2) not in idx
                ):
                    rec[
                        f"{s}_raw_center"
                    ] = np.nan

                    rec[
                        f"{s}_raw_3x3"
                    ] = np.nan

                    rec[
                        f"{s}_balanced_center"
                    ] = np.nan

                    rec[
                        f"{s}_balanced_3x3"
                    ] = np.nan

                    continue

                i = idx[int(rr.start1)]
                j = idx[int(rr.start2)]

                rec[
                    f"{s}_raw_center"
                ] = float(raw[i, j])

                rec[
                    f"{s}_raw_3x3"
                ] = window_sum(
                    raw,
                    i,
                    j,
                    args.quant_radius_bins
                )

                x = bal[i, j]

                rec[
                    f"{s}_balanced_center"
                ] = (
                    float(x)
                    if np.isfinite(x)
                    else np.nan
                )

                rec[
                    f"{s}_balanced_3x3"
                ] = window_sum(
                    bal,
                    i,
                    j,
                    args.quant_radius_bins
                )

            outrows.append(rec)

    out = pd.DataFrame(outrows)

    out.to_csv(
        outdir /
        "pooled_loop_candidates_counts.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip"
    )

    dots.to_csv(
        outdir /
        "pooled_loop_candidates.bedpe.tsv",
        sep="\t",
        index=False
    )

    print(
        "Condition-blind pooled loop candidates: "
        f"{len(out):,}"
    )

    n = (
        2 * args.quant_radius_bins + 1
    )

    print(
        f"Differential quantification uses "
        f"a {n}x{n} window around each "
        f"{args.resolution // 1000}-kb "
        f"loop pixel."
    )


if __name__ == "__main__":
    main()
