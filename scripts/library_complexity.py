#!/usr/bin/env python3

import argparse
import math
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_stats(path):
    stats = {}

    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 2:
                continue

            try:
                stats[parts[0]] = float(parts[1])
            except ValueError:
                pass

    return stats


def get(stats, key):
    if key not in stats:
        raise KeyError(
            f"Missing expected pairtools statistic: {key}"
        )

    return stats[key]


def estimate_complexity(mapped, unique):
    """
    Solve:
        unique = C * (1 - exp(-mapped / C))

    This is the same simple occupancy model underlying
    the pairtools naive-complexity estimate.
    """

    if mapped <= 0 or unique <= 0:
        return np.nan

    if unique >= mapped:
        return np.inf

    def f(c):
        return c * (1.0 - math.exp(-mapped / c)) - unique

    lo = max(unique, 1.0)
    hi = max(mapped, unique) * 10.0

    while f(hi) <= 0 and hi < 1e18:
        hi *= 10.0

    if f(hi) <= 0:
        return np.nan

    for _ in range(200):
        mid = (lo + hi) / 2.0

        if f(mid) > 0:
            hi = mid
        else:
            lo = mid

    return (lo + hi) / 2.0


def projected_unique(mapped_depth, complexity):
    if not np.isfinite(complexity):
        return mapped_depth

    return complexity * (
        1.0 - math.exp(
            -mapped_depth / complexity
        )
    )


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--stats",
        nargs="+",
        required=True
    )

    ap.add_argument(
        "--samples",
        nargs="+",
        required=True
    )

    ap.add_argument(
        "--conditions",
        nargs="+",
        required=True
    )

    ap.add_argument(
        "--multipliers",
        nargs="+",
        type=float,
        required=True
    )

    ap.add_argument(
        "--warning-multiplier",
        type=float,
        required=True
    )

    ap.add_argument(
        "--min-marginal-unique-fraction",
        type=float,
        required=True
    )

    ap.add_argument(
        "--out-tsv",
        required=True
    )

    ap.add_argument(
        "--out-pdf",
        required=True
    )

    ap.add_argument(
        "--out-status",
        required=True
    )

    args = ap.parse_args()

    if not (
        len(args.stats)
        == len(args.samples)
        == len(args.conditions)
    ):
        raise ValueError(
            "stats, samples and conditions "
            "must have equal lengths"
        )

    multipliers = sorted(
        set(
            [1.0]
            + list(args.multipliers)
            + [args.warning_multiplier]
        )
    )

    rows = []

    for path, sample, condition in zip(
        args.stats,
        args.samples,
        args.conditions
    ):
        s = parse_stats(path)

        total = get(s, "total")
        mapped = get(s, "total_mapped")
        dups = get(s, "total_dups")
        nodup = get(s, "total_nodups")

        cis = get(s, "cis")
        trans = get(s, "trans")
        cis1 = get(s, "cis_1kb+")

        valid = cis1 + trans

        complexity = s.get(
            "summary/complexity_naive",
            np.nan
        )

        complexity_source = (
            "pairtools_complexity_naive"
        )

        if (
            not np.isfinite(complexity)
            or complexity <= nodup
        ):
            complexity = estimate_complexity(
                mapped,
                nodup
            )

            complexity_source = (
                "occupancy_model_fallback"
            )

        mapping_fraction = (
            mapped / total
            if total
            else np.nan
        )

        cis_fraction = (
            cis / nodup
            if nodup
            else np.nan
        )

        trans_fraction = (
            trans / nodup
            if nodup
            else np.nan
        )

        valid_fraction = (
            valid / nodup
            if nodup
            else np.nan
        )

        complexity_utilized = (
            nodup / complexity
            if np.isfinite(complexity)
            and complexity > 0
            else np.nan
        )

        for mult in multipliers:
            projected_total = total * mult
            projected_mapped = mapped * mult

            projected_nodup = projected_unique(
                projected_mapped,
                complexity
            )

            projected_nodup = min(
                projected_nodup,
                projected_mapped
            )

            projected_dups = max(
                projected_mapped
                - projected_nodup,
                0
            )

            extra_mapped = max(
                projected_mapped - mapped,
                0
            )

            extra_nodup = max(
                projected_nodup - nodup,
                0
            )

            marginal_unique_fraction = (
                extra_nodup / extra_mapped
                if extra_mapped > 0
                else np.nan
            )

            gain_pct = (
                100.0 * extra_nodup / nodup
                if nodup
                else np.nan
            )

            projected_dup_pct = (
                100.0
                * projected_dups
                / projected_mapped
                if projected_mapped
                else np.nan
            )

            low_return = (
                bool(
                    marginal_unique_fraction
                    < args.min_marginal_unique_fraction
                )
                if np.isfinite(
                    marginal_unique_fraction
                )
                else False
            )

            rows.append({
                "sample": sample,
                "condition": condition,
                "depth_multiplier": mult,

                "observed_total_pairs": total,
                "observed_mapped_pairs": mapped,
                "observed_duplicate_pairs": dups,
                "observed_nodup_pairs": nodup,
                "observed_cis_pairs": cis,
                "observed_trans_pairs": trans,
                "observed_valid_pairs": valid,

                "observed_mapping_fraction":
                    mapping_fraction,

                "observed_duplicate_fraction_mapped":
                    dups / mapped
                    if mapped
                    else np.nan,

                "complexity_naive":
                    complexity,

                "complexity_source":
                    complexity_source,

                "complexity_utilized_fraction":
                    complexity_utilized,

                "projected_total_pairs":
                    projected_total,

                "projected_mapped_pairs":
                    projected_mapped,

                "projected_duplicate_pairs":
                    projected_dups,

                "projected_nodup_pairs":
                    projected_nodup,

                "projected_cis_pairs":
                    projected_nodup
                    * cis_fraction,

                "projected_trans_pairs":
                    projected_nodup
                    * trans_fraction,

                "projected_valid_pairs":
                    projected_nodup
                    * valid_fraction,

                "projected_duplicate_pct_mapped":
                    projected_dup_pct,

                "additional_nodup_pairs":
                    extra_nodup,

                "additional_nodup_gain_pct":
                    gain_pct,

                "marginal_unique_fraction":
                    marginal_unique_fraction,

                "low_return":
                    low_return,
            })

    df = pd.DataFrame(rows)

    os.makedirs(
        os.path.dirname(args.out_tsv),
        exist_ok=True
    )

    df.to_csv(
        args.out_tsv,
        sep="\t",
        index=False
    )

    # --------------------------------------------------------
    # Status at the configured warning multiplier
    # --------------------------------------------------------

    check = df[
        np.isclose(
            df["depth_multiplier"],
            args.warning_multiplier
        )
    ].copy()

    with open(args.out_status, "w") as out:
        out.write(
            "# Naive same-library complexity projection\n"
        )

        out.write(
            "# LOW_RETURN threshold is a configurable QC "
            "heuristic, not a biological cutoff.\n"
        )

        out.write(
            "# Projection assumes unchanged mapping and "
            "contact-composition behaviour.\n\n"
        )

        for _, r in check.iterrows():
            status = (
                "LOW_RETURN"
                if r["low_return"]
                else "SUBSTANTIAL_RETURN"
            )

            out.write(
                f"{r['sample']}\t"
                f"{status}\t"
                f"depth={r['depth_multiplier']:.2f}x\t"
                f"complexity={r['complexity_naive']:.0f}\t"
                f"additional_no_dup="
                f"{r['additional_nodup_pairs']:.0f}\t"
                f"gain_vs_current="
                f"{r['additional_nodup_gain_pct']:.1f}%\t"
                f"marginal_unique="
                f"{100*r['marginal_unique_fraction']:.1f}%\n"
            )

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(12, 4.8)
    )

    for sample, d in df.groupby(
        "sample",
        sort=False
    ):
        axes[0].plot(
            d["projected_mapped_pairs"] / 1e6,
            d["projected_nodup_pairs"] / 1e6,
            marker="o",
            label=sample
        )

        dd = d[
            d["depth_multiplier"] > 1
        ]

        axes[1].plot(
            dd["depth_multiplier"],
            100
            * dd["marginal_unique_fraction"],
            marker="o",
            label=sample
        )

    axes[0].set_xlabel(
        "Mapped read pairs (millions)"
    )

    axes[0].set_ylabel(
        "Expected unique/no-dup pairs (millions)"
    )

    axes[0].set_title(
        "Naive library-complexity projection"
    )

    axes[1].axhline(
        100
        * args.min_marginal_unique_fraction,
        linestyle="--",
        linewidth=1
    )

    axes[1].set_xlabel(
        "Sequencing-depth multiplier"
    )

    axes[1].set_ylabel(
        "Expected unique fraction of additional mapped pairs (%)"
    )

    axes[1].set_title(
        "Marginal yield from additional sequencing"
    )

    for ax in axes:
        ax.spines[
            ["top", "right"]
        ].set_visible(False)

    axes[0].legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(args.out_pdf)
    plt.close(fig)


if __name__ == "__main__":
    main()
