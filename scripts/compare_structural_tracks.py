#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np


def read_e1(path, sample):
    df = pd.read_csv(path, sep="\t")
    needed = {"chrom", "start", "end", "E1"}

    if not needed.issubset(df.columns):
        raise ValueError(f"{path}: missing E1 columns")

    return (
        df[["chrom", "start", "end", "E1"]]
        .rename(columns={"E1": sample})
    )


def read_ins(path, sample, window):
    df = pd.read_csv(path, sep="\t")

    preferred = f"log2_insulation_score_{window}"

    if preferred in df.columns:
        c = preferred
    else:
        candidates = [
            x for x in df.columns
            if x.startswith("log2_insulation_score")
            and str(window) in x
        ]

        if len(candidates) != 1:
            raise ValueError(
                f"{path}: cannot uniquely identify "
                f"log2 insulation score for window {window}. "
                f"Candidates: {candidates}"
            )

        c = candidates[0]

    return (
        df[["chrom", "start", "end", c]]
        .rename(columns={c: sample})
    )


def merge_tracks(paths, samples, reader):
    if len(paths) != len(samples):
        raise ValueError(
            "Number of track files does not match number of samples"
        )

    out = None

    for path, sample in zip(paths, samples):
        x = reader(path, sample)

        if out is None:
            out = x
        else:
            out = out.merge(
                x,
                on=["chrom", "start", "end"],
                how="outer"
            )

    return out


def add_effects(df, control_samples, test_samples):
    missing = [
        s for s in control_samples + test_samples
        if s not in df.columns
    ]

    if missing:
        raise ValueError(
            "Samples missing from merged track table: "
            + ", ".join(missing)
        )

    # Group means
    df["control_mean"] = df[control_samples].mean(
        axis=1,
        skipna=True
    )

    df["test_mean"] = df[test_samples].mean(
        axis=1,
        skipna=True
    )

    df["test_minus_control_mean"] = (
        df["test_mean"] - df["control_mean"]
    )

    # Per-test-sample effects relative to control mean
    delta_cols = []

    for sample in test_samples:
        col = f"{sample}_minus_control_mean"

        df[col] = (
            df[sample] - df["control_mean"]
        )

        delta_cols.append(col)

    # Require all test samples to have non-missing effects
    complete = df[delta_cols].notna().all(axis=1)

    # Same non-zero sign in every test sample
    signs = np.sign(df[delta_cols])

    same_positive = (signs > 0).all(axis=1)
    same_negative = (signs < 0).all(axis=1)

    df["all_test_same_direction"] = (
        complete &
        (same_positive | same_negative)
    )

    df["consensus_direction"] = np.select(
        [
            complete & same_positive,
            complete & same_negative
        ],
        [
            "up",
            "down"
        ],
        default="discordant"
    )

    return df


def write_table(df, path):
    compression = (
        "gzip"
        if str(path).endswith(".gz")
        else None
    )

    df.to_csv(
        path,
        sep="\t",
        index=False,
        compression=compression
    )


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--e1",
        nargs="+",
        required=True
    )

    ap.add_argument(
        "--insulation",
        nargs="+",
        required=True
    )

    ap.add_argument(
        "--samples",
        nargs="+",
        required=True
    )

    ap.add_argument(
        "--control-samples",
        nargs="+",
        required=True
    )

    ap.add_argument(
        "--test-samples",
        nargs="+",
        required=True
    )

    ap.add_argument(
        "--insulation-window",
        type=int,
        required=True
    )

    ap.add_argument(
        "--out-compartments",
        required=True
    )

    ap.add_argument(
        "--out-insulation",
        required=True
    )

    args = ap.parse_args()

    if not (
        len(args.e1)
        == len(args.insulation)
        == len(args.samples)
    ):
        raise ValueError(
            "Track lists and samples must have equal lengths"
        )

    if len(args.control_samples) < 2:
        raise ValueError(
            "At least two control samples are required"
        )

    if len(args.test_samples) < 2:
        raise ValueError(
            "At least two test samples are required"
        )

    unknown = (
        set(args.control_samples + args.test_samples)
        - set(args.samples)
    )

    if unknown:
        raise ValueError(
            "Control/test samples absent from --samples: "
            + ", ".join(sorted(unknown))
        )

    comp = merge_tracks(
        args.e1,
        args.samples,
        read_e1
    )

    comp = add_effects(
        comp,
        args.control_samples,
        args.test_samples
    )

    ins = merge_tracks(
        args.insulation,
        args.samples,
        lambda path, sample: read_ins(
            path,
            sample,
            args.insulation_window
        )
    )

    ins = add_effects(
        ins,
        args.control_samples,
        args.test_samples
    )

    write_table(
        comp,
        args.out_compartments
    )

    write_table(
        ins,
        args.out_insulation
    )


if __name__ == "__main__":
    main()
