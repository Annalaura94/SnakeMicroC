#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import cooler
import numpy as np
import pandas as pd


def parse_csv(x):
    return [
        v for v in x.split(",")
        if v
    ]


def window_sum(mat, i, j, r):
    i0 = max(0, i - r)
    i1 = min(mat.shape[0], i + r + 1)

    j0 = max(0, j - r)
    j1 = min(mat.shape[1], j + r + 1)

    return float(
        mat[
            i0:i1,
            j0:j1
        ].sum()
    )


def finite_or_nan(x):
    x = float(x)

    if np.isfinite(x):
        return x

    return np.nan


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--samples",
        required=True
    )

    ap.add_argument(
        "--resolution",
        type=int,
        required=True
    )

    ap.add_argument(
        "--quant-radius-bins",
        type=int,
        required=True
    )

    ap.add_argument(
        "--loops",
        required=True
    )

    ap.add_argument(
        "--cool-pattern",
        required=True
    )

    ap.add_argument(
        "--outdir",
        required=True
    )

    args = ap.parse_args()

    samples = parse_csv(
        args.samples
    )

    if len(samples) < 2:
        raise ValueError(
            "At least two samples are required"
        )

    outdir = Path(
        args.outdir
    )

    outdir.mkdir(
        parents=True,
        exist_ok=True
    )

    loops = pd.read_csv(
        args.loops,
        sep="\t"
    )

    required = [
        "loop_id",
        "chrom1",
        "start1",
        "end1",
        "chrom2",
        "start2",
        "end2",
        "distance_bp",
    ]

    missing = [
        x for x in required
        if x not in loops.columns
    ]

    if missing:
        raise ValueError(
            f"Loop table missing {missing}; "
            f"columns={list(loops.columns)}"
        )

    if loops["loop_id"].duplicated().any():
        raise ValueError(
            "loop_id values are not unique"
        )

    # --------------------------------------------------------
    # Open sample coolers and record total contact depth.
    # These totals are used as edgeR exposure/library sizes.
    # --------------------------------------------------------

    cools = {}

    library_rows = []

    for sample in samples:
        path = args.cool_pattern.format(
            sample=sample,
            resolution=args.resolution
        )

        clr = cooler.Cooler(
            path
        )

        if clr.binsize != args.resolution:
            raise ValueError(
                f"{sample}: cooler binsize "
                f"{clr.binsize} != expected "
                f"{args.resolution}"
            )

        total = clr.info.get(
            "sum",
            None
        )

        if total is None:
            raise ValueError(
                f"{sample}: cooler info does not "
                "contain total contact sum"
            )

        total = int(total)

        if total <= 0:
            raise ValueError(
                f"{sample}: invalid cooler sum {total}"
            )

        cools[sample] = clr

        library_rows.append({
            "sample": sample,
            "total_contacts": total,
            "resolution": args.resolution,
            "cooler": path,
        })

    library_sizes = pd.DataFrame(
        library_rows
    )

    library_sizes.to_csv(
        outdir /
        "loop_library_sizes.tsv",
        sep="\t",
        index=False
    )

    # --------------------------------------------------------
    # Quantify every condition-blind candidate in every sample
    # --------------------------------------------------------

    outrows = []

    for chrom, d in loops.groupby(
        "chrom1",
        sort=False
    ):
        cache = {}

        for sample, clr in cools.items():

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

            cache[sample] = (
                start_to_idx,
                raw,
                bal
            )

        for rr in d.itertuples(
            index=False
        ):
            rec = rr._asdict()

            for sample in samples:
                idx, raw, bal = cache[
                    sample
                ]

                start1 = int(
                    rr.start1
                )

                start2 = int(
                    rr.start2
                )

                if (
                    start1 not in idx
                    or
                    start2 not in idx
                ):
                    raise ValueError(
                        f"{rr.loop_id}: anchor missing "
                        f"from {sample} cooler: "
                        f"{chrom}:{start1}, "
                        f"{chrom}:{start2}"
                    )

                i = idx[start1]
                j = idx[start2]

                raw_center = float(
                    raw[i, j]
                )

                raw_count = window_sum(
                    raw,
                    i,
                    j,
                    args.quant_radius_bins
                )

                bal_center = finite_or_nan(
                    bal[i, j]
                )

                bal_count = finite_or_nan(
                    window_sum(
                        bal,
                        i,
                        j,
                        args.quant_radius_bins
                    )
                )

                # Raw contact counts must remain integer.
                rec[
                    f"{sample}_raw_center"
                ] = int(
                    round(raw_center)
                )

                rec[
                    f"{sample}_raw_count"
                ] = int(
                    round(raw_count)
                )

                # Balanced values are QC/visualization only.
                rec[
                    f"{sample}_balanced_center"
                ] = bal_center

                rec[
                    f"{sample}_balanced_count"
                ] = bal_count

            outrows.append(rec)

    out = pd.DataFrame(
        outrows
    )

    out.to_csv(
        outdir /
        "pooled_loop_candidates_counts.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip"
    )

    n = (
        2 * args.quant_radius_bins + 1
    )

    print(
        "Condition-blind pooled loop candidates: "
        f"{len(out):,}"
    )

    print(
        "Differential quantification uses raw "
        f"{n}x{n} contact windows around each "
        f"{args.resolution // 1000}-kb loop pixel."
    )

    print(
        "Library sizes for edgeR are total "
        "contacts stored in the corresponding "
        "full-depth sample coolers."
    )


if __name__ == "__main__":
    main()
