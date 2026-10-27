#!/usr/bin/env python3

from __future__ import annotations

import argparse
import html
from pathlib import Path

import numpy as np
import pandas as pd


def read_tsv(path):
    return pd.read_csv(path, sep="\t")


def metric_dict(df):
    if not {"metric", "value"}.issubset(df.columns):
        raise ValueError("Expected columns: metric, value")
    return dict(zip(df["metric"], df["value"]))


def fmt_int(x):
    return f"{int(x):,}"


def fmt_float(x, n=3):
    return f"{float(x):.{n}f}"


def table_md(df):
    if df.empty:
        return "_No entries._"
    return df.to_markdown(index=False)


def interpretation(
    candidates,
    tested,
    fdr_sig,
    effect_sig,
    core,
    increased,
    decreased,
    norm,
    fdr_cut,
    test_label,
):
    lines = []

    retained_pct = (
        100.0 * tested / candidates
        if candidates > 0 else np.nan
    )

    lines.append(
        f"Mustache identified **{fmt_int(candidates)} "
        f"condition-blind candidate loops**. "
        f"After independent low-count filtering, "
        f"**{fmt_int(tested)} ({retained_pct:.1f}%)** "
        f"were retained for differential testing."
    )

    if fdr_sig == 0:
        lines.append(
            f"No loop reached FDR < {fdr_cut:g}. "
            "At this resolution and sequencing depth, "
            "there is therefore no statistically supported "
            "genome-wide differential-loop signal."
        )
    else:
        lines.append(
            f"**{fmt_int(fdr_sig)} loops** reached FDR < {fdr_cut:g}."
        )

    if effect_sig == 0:
        lines.append(
            "None of the tested loops simultaneously passed "
            "the pre-specified FDR and absolute log2FC threshold."
        )
    else:
        lines.append(
            f"Applying the pre-specified effect-size criterion "
            f"identified **{fmt_int(effect_sig)} differential loops**."
        )

    if core == 0:
        lines.append(
            "No differential loop additionally showed the required "
            "directional concordance across all test samples. "
            "Thus, no loop satisfies the stringent concordant-core "
            "definition."
        )
    else:
        lines.append(
            f"**{fmt_int(core)} loops** also changed in the same "
            f"direction across all {test_label} samples and constitute the "
            f"concordant core: **{fmt_int(increased)} increased** "
            f"and **{fmt_int(decreased)} decreased** in {test_label}."
        )

    if not norm.empty:
        mn = norm["TMM_norm_factor"].min()
        mx = norm["TMM_norm_factor"].max()

        if mn >= 0.8 and mx <= 1.25:
            lines.append(
                f"TMM normalization factors ranged from "
                f"{mn:.3f} to {mx:.3f}, indicating relatively "
                f"modest composition scaling among samples."
            )
        else:
            lines.append(
                f"TMM normalization factors ranged from "
                f"{mn:.3f} to {mx:.3f}. "
                f"The wider scaling range should be considered "
                f"when interpreting global differences in loop signal."
            )

    return "\n\n".join(lines)


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--summary", required=True)
    ap.add_argument("--filtering", required=True)
    ap.add_argument("--normalization", required=True)
    ap.add_argument("--all-results", required=True)
    ap.add_argument("--core-results", required=True)
    ap.add_argument("--library-sizes", required=True)
    ap.add_argument("--parameters", required=True)
    ap.add_argument("--versions", required=True)
    ap.add_argument("--resolution", type=int, required=True)
    ap.add_argument("--fdr", type=float, required=True)
    ap.add_argument("--lfc", type=float, required=True)
    ap.add_argument("--output-md", required=True)
    ap.add_argument("--output-html", required=True)
    ap.add_argument("--quant-radius", type=int, default=1)
    ap.add_argument("--control-label", default="control")
    ap.add_argument("--test-label", default="test")

    args = ap.parse_args()

    summary = read_tsv(args.summary)
    filtering = read_tsv(args.filtering)
    norm = read_tsv(args.normalization)
    res = read_tsv(args.all_results)
    core = read_tsv(args.core_results)
    libs = read_tsv(args.library_sizes)
    params = read_tsv(args.parameters)

    versions_text = Path(args.versions).read_text().strip()

    sm = metric_dict(summary)
    fm = metric_dict(filtering)

    candidates = int(
        sm.get(
            "condition_blind_candidates",
            fm.get("condition_blind_candidates", len(res))
        )
    )

    tested = int(
        sm.get("tested_after_filter", len(res))
    )

    fdr_sig = int(
        sm.get("FDR_significant", 0)
    )

    effect_sig = int(
        sm.get("FDR_and_effect_size", 0)
    )

    n_core = int(
        sm.get("concordant_core", len(core))
    )

    increased = int(
        sm.get("core_increased", 0)
    )

    decreased = int(
        sm.get("core_decreased", 0)
    )

    interpretation_text = interpretation(
        candidates,
        tested,
        fdr_sig,
        effect_sig,
        n_core,
        increased,
        decreased,
        norm,
        args.fdr,
        args.test_label,
    )

    # --------------------------------------------------------
    # Top differential loops
    # --------------------------------------------------------

    top_cols = [
        x for x in [
            "loop_id",
            "chrom1",
            "start1",
            "chrom2",
            "start2",
            "distance_bp",
            "logFC",
            "logCPM",
            "PValue",
            "FDR",
            "differential_loop",
            "concordant_core"
        ]
        if x in res.columns
    ]

    top = (
        res
        .sort_values(
            ["FDR", "PValue"],
            na_position="last"
        )
        .head(20)[top_cols]
        .copy()
    )

    for col in [
        "logFC",
        "logCPM",
        "PValue",
        "FDR"
    ]:
        if col in top.columns:
            top[col] = top[col].map(
                lambda x:
                f"{x:.4g}"
                if pd.notna(x)
                else ""
            )

    # --------------------------------------------------------
    # Library table
    # --------------------------------------------------------

    libs_display = libs.copy()

    if "total_contacts" in libs_display.columns:
        libs_display["total_contacts"] = (
            libs_display["total_contacts"]
            .map(lambda x: f"{int(x):,}")
        )

    norm_display = norm.copy()

    for col in [
        "library_size",
        "effective_library_size"
    ]:
        if col in norm_display.columns:
            norm_display[col] = (
                norm_display[col]
                .map(lambda x: f"{float(x):,.0f}")
            )

    if "TMM_norm_factor" in norm_display.columns:
        norm_display["TMM_norm_factor"] = (
            norm_display["TMM_norm_factor"]
            .map(lambda x: f"{float(x):.4f}")
        )

    params_display = params.copy()

    # --------------------------------------------------------
    # Methods paragraph
    # --------------------------------------------------------

    kb = args.resolution // 1000
    window = 2 * args.quant_radius + 1
    condition_comparison = f"{args.test_label} versus {args.control_label}"

    methods = (
        f"Loop candidates were identified condition-blind from the "
        f"pooled Micro-C contact map using Mustache at {kb}-kb "
        f"resolution. Candidate loops were restricted to the "
        f"pre-specified genomic distance interval and quantified "
        f"independently in each biological sample using raw contact "
        f"counts from a {window}x{window}-bin window centred on each candidate. "
        f"Differential loop strength between {args.control_label} and {args.test_label} samples "
        f"was tested with edgeR quasi-likelihood negative-binomial "
        f"models. Full-depth contact-map counts were supplied as "
        f"library-size exposure values and TMM factors were used for "
        f"between-sample composition normalization. Differential "
        f"loops were pre-defined as FDR < {args.fdr:g} and "
        f"|log2FC| >= {args.lfc:g}. A stringent concordant core "
        f"additionally required all test samples to change in the same "
        f"direction as the group-level effect."
    )

    # --------------------------------------------------------
    # Markdown
    # --------------------------------------------------------

    md = f"""# Micro-C differential loop analysis

## Analysis overview

**Resolution:** {kb} kb

**Condition comparison:** {condition_comparison}

**Candidate caller:** Mustache

**Candidate loops:** {fmt_int(candidates)}

**Loops tested after filtering:** {fmt_int(tested)}

**FDR-significant loops:** {fmt_int(fdr_sig)}

**FDR + effect-size loops:** {fmt_int(effect_sig)}

**Concordant core loops:** {fmt_int(n_core)}

---

## 1. What was tested and why

Loop discovery was performed **condition-blind**, so condition labels were not used to define the candidate-loop universe. This separates feature discovery from the subsequent differential test.

Candidate loops were then quantified separately in each biological sample. The inferential comparison therefore retains the biological replication structure rather than comparing pooled WT and pooled KO maps.

## 2. Mustache parameters

{table_md(params_display)}

## 3. Sequencing/contact depth used for quantification

{table_md(libs_display)}

The differential model uses full-depth sample contact counts as library-size exposure values. Balanced contact values are retained for visualization/QC but are not used as edgeR counts.

## 4. Filtering and normalization

{table_md(filtering)}

### Normalization factors

{table_md(norm_display)}

Low-count loops are removed independently of differential status before model fitting. TMM normalization then accounts for remaining composition differences among samples.

## 5. Differential results

{table_md(summary)}

### Interpretation

{interpretation_text}

## 6. Strongest statistical signals

The table below shows the 20 loops with the smallest FDR values, independently of whether they pass the final pre-specified thresholds.

{table_md(top)}

## 7. Interpretation limits

A statistically differential contact does not by itself identify the regulatory target or demonstrate enhancer-promoter causality. Biological interpretation should therefore be performed only after genomic annotation of loop anchors and integration with gene expression and regulatory-element information.

Per-test-sample fold changes are used as a **descriptive concordance criterion**, not as independent n=1 statistical tests.

## 8. Methods-ready description

{methods}

## 9. Software provenance

```text
{versions_text}
```

## 10. Associated QC

The analysis directory also contains:

- `EXP_differential_loops_volcano.pdf`
- `EXP_differential_loops_QC.pdf`
- `loop_normalization_factors.tsv`
- `loop_filtering_summary.tsv`
- `differential_loops_all.tsv.gz`
- `differential_loops_FDR_effect.tsv`
- `differential_loops_core.tsv`
- `mustache_environment_versions.txt`

These files should be retained together with this report as the provenance package for the analysis.
"""

    Path(args.output_md).write_text(md)

    # --------------------------------------------------------
    # HTML
    # --------------------------------------------------------

    def html_table(df):
        if df.empty:
            return "<p><em>No entries.</em></p>"

        return df.to_html(
            index=False,
            border=0,
            classes="dataframe",
            escape=True
        )

    interp_html = "".join(
        f"<p>{html.escape(x)}</p>"
        for x in interpretation_text.split("\n\n")
    )

    html_doc = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>Micro-C differential loop analysis</title>
<style>
body {{
    font-family: Arial, Helvetica, sans-serif;
    max-width: 1150px;
    margin: 40px auto;
    padding: 0 30px 60px 30px;
    line-height: 1.55;
    color: #222;
}}
h1 {{
    border-bottom: 3px solid #222;
    padding-bottom: 12px;
}}
h2 {{
    margin-top: 38px;
    border-bottom: 1px solid #bbb;
    padding-bottom: 6px;
}}
.summary {{
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 12px;
    margin: 24px 0;
}}
.card {{
    border: 1px solid #ccc;
    border-radius: 7px;
    padding: 14px;
}}
.card .number {{
    font-size: 1.7em;
    font-weight: bold;
}}
table.dataframe {{
    border-collapse: collapse;
    width: 100%;
    font-size: 0.90em;
    margin: 15px 0 25px 0;
}}
table.dataframe th,
table.dataframe td {{
    border-bottom: 1px solid #ddd;
    padding: 7px 9px;
    text-align: left;
}}
table.dataframe th {{
    background: #f3f3f3;
}}
.note {{
    background: #f7f7f7;
    border-left: 4px solid #555;
    padding: 12px 16px;
    margin: 18px 0;
}}
code {{
    background: #f1f1f1;
    padding: 2px 4px;
}}
</style>
</head>
<body>

<h1>Micro-C differential loop analysis</h1>

<p><strong>Resolution:</strong> {kb} kb<br>
<strong>Comparison:</strong> {condition_comparison}<br>
<strong>Candidate caller:</strong> Mustache</p>

<div class="summary">
<div class="card"><div class="number">{fmt_int(candidates)}</div>candidate loops</div>
<div class="card"><div class="number">{fmt_int(tested)}</div>loops tested</div>
<div class="card"><div class="number">{fmt_int(fdr_sig)}</div>FDR significant</div>
<div class="card"><div class="number">{fmt_int(effect_sig)}</div>FDR + effect size</div>
<div class="card"><div class="number">{fmt_int(n_core)}</div>concordant core</div>
<div class="card"><div class="number">{fmt_int(increased)} / {fmt_int(decreased)}</div>core increased / decreased</div>
</div>

<h2>1. Analysis design</h2>

<p>Loop discovery was performed <strong>condition-blind</strong>. Candidate loops were subsequently quantified independently in each biological sample, preserving biological replication for the differential analysis.</p>

<h2>2. Mustache parameters</h2>
{html_table(params_display)}

<h2>3. Contact depth</h2>
{html_table(libs_display)}

<p>Raw contact counts are used for differential inference. Balanced values are retained only for QC and visualization.</p>

<h2>4. Filtering and normalization</h2>
{html_table(filtering)}

<h3>Normalization factors</h3>
{html_table(norm_display)}

<h2>5. Differential-loop results</h2>
{html_table(summary)}

<div class="note">
{interp_html}
</div>

<p>
<a href="EXP_differential_loops_volcano.pdf">Open volcano plot</a> |
<a href="EXP_differential_loops_QC.pdf">Open differential-loop QC</a>
</p>

<h2>6. Strongest statistical signals</h2>
<p>Top 20 loops ranked by FDR. Ranking does not replace the pre-specified significance criteria.</p>
{html_table(top)}

<h2>7. Interpretation limits</h2>

<p>A differential contact does not by itself establish enhancer-promoter causality or identify the regulated gene. Anchor annotation and integration with expression/regulatory datasets are required for biological interpretation.</p>

<p>Per-test-sample effects are used for descriptive directional concordance, not for independent per-sample differential tests.</p>

<h2>8. Methods-ready description</h2>

<div class="note">
<p>{html.escape(methods)}</p>
</div>

<h2>9. Software provenance</h2>

<pre>{html.escape(versions_text)}</pre>

<h2>10. Provenance package</h2>

<p>The complete analysis directory contains the candidate-loop set, raw per-sample quantification, filtering summary, normalization factors, full edgeR results, concordant-core results, QC figures and software-version information.</p>

</body>
</html>
"""

    Path(args.output_html).write_text(
        html_doc
    )

    print(
        f"Report written: {args.output_html}"
    )
    print(
        f"Markdown written: {args.output_md}"
    )


if __name__ == "__main__":
    main()
