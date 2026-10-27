#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--view", required=True)

    ap.add_argument("--resolution", type=int, required=True)
    ap.add_argument("--min-distance-bp", type=int, required=True)
    ap.add_argument("--max-distance-bp", type=int, required=True)

    args = ap.parse_args()

    df = pd.read_csv(
        args.input,
        sep="\t"
    )

    df.columns = [
        str(x).strip().lstrip("#")
        for x in df.columns
    ]

    aliases = {
        "BIN1_CHR": "chrom1",
        "BIN1_START": "start1",
        "BIN1_END": "end1",
        "BIN2_CHROMOSOME": "chrom2",
        "BIN2_CHR": "chrom2",
        "BIN2_START": "start2",
        "BIN2_END": "end2",
        "FDR": "mustache_fdr",
        "DETECTION_SCALE": "mustache_detection_scale",
    }

    df = df.rename(
        columns={
            k: v
            for k, v in aliases.items()
            if k in df.columns
        }
    )

    required = [
        "chrom1",
        "start1",
        "end1",
        "chrom2",
        "start2",
        "end2",
    ]

    missing = [
        x for x in required
        if x not in df.columns
    ]

    if missing:
        raise ValueError(
            f"Missing Mustache columns after renaming: {missing}. "
            f"Observed columns: {list(df.columns)}"
        )

    coord_cols = [
        "start1",
        "end1",
        "start2",
        "end2",
    ]

    for col in coord_cols:
        df[col] = pd.to_numeric(
            df[col],
            errors="coerce"
        )

    df = df.dropna(
        subset=coord_cols
    ).copy()

    for col in coord_cols:
        df[col] = df[col].astype(np.int64)

    df["chrom1"] = df["chrom1"].astype(str)
    df["chrom2"] = df["chrom2"].astype(str)

    # --------------------------------------------------------
    # Restrict to chromosomes/regions defined by MAIN_VIEW
    # --------------------------------------------------------

    view = pd.read_csv(
        args.view,
        sep="\t",
        comment="#",
        header=None,
        dtype=str
    )

    chrom_order = []

    for chrom in view.iloc[:, 0].tolist():
        if chrom.lower() in {
            "chrom",
            "chromosome"
        }:
            continue

        if chrom not in chrom_order:
            chrom_order.append(chrom)

    allowed = set(chrom_order)

    df = df[
        df["chrom1"].isin(allowed) &
        df["chrom2"].isin(allowed)
    ].copy()

    # --------------------------------------------------------
    # Cis loops only
    # --------------------------------------------------------

    df = df[
        df["chrom1"] == df["chrom2"]
    ].copy()

    # --------------------------------------------------------
    # Put upstream anchor first
    # --------------------------------------------------------

    swap = df["start1"] > df["start2"]

    if swap.any():
        for a, b in [
            ("chrom1", "chrom2"),
            ("start1", "start2"),
            ("end1", "end2"),
        ]:
            tmp = df.loc[swap, a].copy()
            df.loc[swap, a] = df.loc[swap, b].values
            df.loc[swap, b] = tmp.values

    df["distance_bp"] = (
        df["start2"] -
        df["start1"]
    ).abs()

    df = df[
        df["distance_bp"].between(
            args.min_distance_bp,
            args.max_distance_bp
        )
    ].copy()

    # --------------------------------------------------------
    # Basic coordinate validation
    # --------------------------------------------------------

    df = df[
        (df["start1"] >= 0) &
        (df["start2"] >= 0) &
        (df["end1"] > df["start1"]) &
        (df["end2"] > df["start2"])
    ].copy()

    if "mustache_fdr" in df.columns:
        df["mustache_fdr"] = pd.to_numeric(
            df["mustache_fdr"],
            errors="coerce"
        )

    if "mustache_detection_scale" in df.columns:
        df["mustache_detection_scale"] = pd.to_numeric(
            df["mustache_detection_scale"],
            errors="coerce"
        )

    # --------------------------------------------------------
    # Exact duplicate removal.
    # If duplicated coordinates exist, retain lowest FDR.
    # --------------------------------------------------------

    coord_key = [
        "chrom1",
        "start1",
        "end1",
        "chrom2",
        "start2",
        "end2",
    ]

    sort_cols = coord_key.copy()

    if "mustache_fdr" in df.columns:
        sort_cols.append("mustache_fdr")

    df = (
        df
        .sort_values(sort_cols)
        .drop_duplicates(
            subset=coord_key,
            keep="first"
        )
        .copy()
    )

    order_map = {
        chrom: i
        for i, chrom in enumerate(chrom_order)
    }

    df["_chrom_order"] = (
        df["chrom1"]
        .map(order_map)
        .fillna(10**9)
    )

    df = (
        df
        .sort_values([
            "_chrom_order",
            "start1",
            "start2",
        ])
        .drop(
            columns="_chrom_order"
        )
        .reset_index(drop=True)
    )

    df.insert(
        0,
        "loop_id",
        [
            f"LOOP_{i:06d}"
            for i in range(
                1,
                len(df) + 1
            )
        ]
    )

    preferred = [
        "loop_id",
        "chrom1",
        "start1",
        "end1",
        "chrom2",
        "start2",
        "end2",
        "distance_bp",
        "mustache_fdr",
        "mustache_detection_scale",
    ]

    preferred = [
        x for x in preferred
        if x in df.columns
    ]

    extra = [
        x for x in df.columns
        if x not in preferred
    ]

    df = df[
        preferred + extra
    ]

    Path(args.output).parent.mkdir(
        parents=True,
        exist_ok=True
    )

    df.to_csv(
        args.output,
        sep="\t",
        index=False
    )

    print(
        f"Canonical Mustache candidates: {len(df):,}"
    )

    print(
        f"Resolution: {args.resolution:,} bp"
    )

    print(
        "Distance range: "
        f"{args.min_distance_bp:,}-"
        f"{args.max_distance_bp:,} bp"
    )


if __name__ == "__main__":
    main()
