#!/usr/bin/env python3

from __future__ import annotations

import argparse
import gzip
import hashlib
import re
import subprocess
from pathlib import Path

import yaml


def first_read_length(path):
    path = str(path)

    # FASTQs may be archived after a completed mapping run. Method generation
    # must not fail or invent read length when those paths are unavailable.
    if not Path(path).is_file():
        return None
    opener = gzip.open if path.endswith(".gz") else open

    with opener(path, "rt") as fh:
        header = fh.readline()
        seq = fh.readline().strip()

    if not header.startswith("@") or not seq:
        return None

    return len(seq)


def fmt_bp(x):
    x = int(x)

    if x % 1_000_000 == 0:
        return f"{x // 1_000_000} Mb"

    if x % 1_000 == 0:
        return f"{x // 1_000} kb"

    return f"{x} bp"


def join_bp(xs):
    vals = [fmt_bp(x) for x in xs]
    if len(vals) <= 1:
        return "".join(vals)
    return ", ".join(vals[:-1]) + " and " + vals[-1]


def get_git_info(pipeline_dir):
    pipeline_dir = str(pipeline_dir)

    try:
        commit = subprocess.check_output(
            [
                "git",
                "-C",
                pipeline_dir,
                "rev-parse",
                "HEAD",
            ],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        commit = "not_available"

    try:
        dirty = subprocess.check_output(
            [
                "git",
                "-C",
                pipeline_dir,
                "status",
                "--porcelain",
            ],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()

        status = (
            "modified working tree"
            if dirty
            else "clean working tree"
        )

    except Exception:
        status = "not_available"

    return commit, status


def sha256(path):
    h = hashlib.sha256()

    with open(path, "rb") as fh:
        for chunk in iter(
            lambda: fh.read(1024 * 1024),
            b"",
        ):
            h.update(chunk)

    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--qc-only", action="store_true",
                    help="Document shallow pilot QC; do not claim deep analyses")
    ap.add_argument("--no-differential", action="store_true",
                    help="Document QC and descriptive analyses only")

    ap.add_argument(
        "--config",
        required=True
    )

    ap.add_argument(
        "--pipeline-dir",
        required=True
    )

    ap.add_argument(
        "--output-md",
        required=True
    )

    ap.add_argument(
        "--output-txt",
        required=True
    )

    args = ap.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    ref = cfg["reference"]
    processing = cfg["processing"]
    qc = cfg["qc"]
    cmaps = cfg["contact_maps"]
    rep = cfg["reproducibility"]
    comparison = cfg.get("comparison") or {}
    samples = cfg["samples"]
    compartments = cfg["compartments"]
    boundaries = cfg["boundaries"]
    contacts = cfg["contacts"]
    loops = cfg["loops"]

    # --------------------------------------------------------
    # Sequencing layout / read length
    # --------------------------------------------------------

    read_lengths = {}

    for sample, dat in samples.items():
        r1_len = first_read_length(dat["r1"])
        r2_len = first_read_length(dat["r2"])

        read_lengths[sample] = (
            r1_len,
            r2_len
        )

    layouts = sorted(
        set(read_lengths.values())
    )

    if (
        len(layouts) == 1
        and layouts[0][0] is not None
        and layouts[0][1] is not None
    ):
        r1_len, r2_len = layouts[0]

        if r1_len == r2_len:
            seq_text = (
                f"paired-end {r1_len}-bp sequencing"
            )
        else:
            seq_text = (
                f"paired-end sequencing "
                f"(R1 {r1_len} bp; R2 {r2_len} bp)"
            )
    else:
        seq_text = "paired-end sequencing"

    # --------------------------------------------------------
    # Experimental design
    # --------------------------------------------------------

    controls = [
        s for s, dat in samples.items()
        if dat["condition"]
        == comparison.get("control_condition")
    ]

    tests = [
        s for s, dat in samples.items()
        if dat["condition"]
        == comparison.get("test_condition")
    ]

    control_names = ", ".join(controls)
    test_names = ", ".join(tests)
    sample_count = len(samples)
    control_count = len(controls)
    test_count = len(tests)
    cond_count = {}
    for dat in samples.values():
        name = dat["condition"]
        cond_count[name] = cond_count.get(name, 0) + 1
    if args.no_differential:
        design_text = (
            f"Micro-C data were analyzed from {sample_count} libraries "
            f"({', '.join(f'{k}: {v}' for k, v in cond_count.items())}). "
            f"Libraries were analyzed using {{seq_text}}. "
            "The configured design does not support an automated group-level "
            "differential analysis, or it was disabled."
        )
    else:
        design_text = None

    # --------------------------------------------------------
    # Loop parameters
    # --------------------------------------------------------

    loop_res = int(
        loops["primary_resolution"]
    )

    radius = int(
        loops["quant_radius_bins"]
    )

    window_n = (
        2 * radius + 1
    )

    sensitivity_res = loops.get(
        "sensitivity_resolutions",
        []
    )

    # --------------------------------------------------------
    # Provenance
    # --------------------------------------------------------

    commit, git_status = get_git_info(
        args.pipeline_dir
    )

    config_hash = sha256(
        args.config
    )

    if args.qc_only:
        groups = ", ".join(f"{k} ({v['condition']})" for k, v in samples.items())
        parts = [
            ("Design and scope", f"Shallow-sequencing quality control was performed on {sample_count} Micro-C libraries: {groups}. Libraries were assessed using {seq_text}. These reads represent a pilot subsample and cannot establish deep-sequencing QC pass/fail, compartment changes, domain boundaries or differential chromatin interactions."),
            ("Mapping and deduplication", f"Paired reads were aligned to {ref['assembly']} using BWA-MEM -5SP -T0, parsed with pairtools (MAPQ >= {processing['mapq']}, walks policy {processing['walks_policy']}, maximum inter-alignment gap {processing['max_inter_align_gap']} bp), sorted, deduplicated, and summarized using pairtools statistics."),
            ("Pilot QC measurements", f"The pilot report summarizes mapping, duplicate rates, cis/trans partition, cis contacts >=1 kb, distance-decay P(s) and library-complexity extrapolation. A cis >=1 kb fraction of {qc['min_cis_ge_1kb_fraction']:.2f} was used as a screening flag, not as a definitive library acceptance threshold. Projected complexity under deeper sequencing is an occupancy-model estimate with substantial uncertainty at shallow depth."),
            ("Exploratory QC maps", f"Raw, unbalanced contact matrices at {fmt_bp(qc['qc_map_resolution'])} resolution were generated directly from pairs and visualized for gross library QC only. Sample correlations were calculated from depth-scaled unbalanced maps. No deep multi-resolution matrix balancing, HiCRep, A/B-compartment calling, insulation analysis, Mustache calling or differential testing was performed in this QC-only run."),
            ("Provenance", f"SnakeMicroC commit {commit} ({git_status}); configuration SHA256 {config_hash}. Software versions are recorded in the report provenance."),
        ]
        md = "# Shallow Micro-C QC - methods and provenance\n\n" + "\n\n".join(f"## {title}\n\n{body}" for title, body in parts) + "\n"
        txt = "Shallow Micro-C QC methods\n\n" + "\n\n".join(body for _, body in parts) + "\n"
        Path(args.output_md).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output_md).write_text(md, encoding="utf8")
        Path(args.output_txt).write_text(txt, encoding="utf8")
        return

    # --------------------------------------------------------
    # Methods text
    # --------------------------------------------------------

    sections = []

    sections.append(
        (
            "Experimental design and sequencing",
            (
                f"Micro-C data were analyzed from {sample_count} libraries: "
                f"{control_count} {comparison.get('control_condition', '')} "
                f"({control_names}) and {test_count} "
                f"{comparison.get('test_condition', '')} "
                f"({test_names}). Libraries were "
                f"analyzed using {seq_text}. The primary comparison "
                f"was {comparison.get('test_condition', '')} versus "
                f"{comparison.get('control_condition', '')}."
            )
        )
    )

    if args.no_differential:
        sections[-1] = ("Experimental design and sequencing", design_text.replace("{seq_text}", seq_text))

    sections.append(
        (
            "Read mapping and contact-pair processing",
            (
                f"Reads were aligned to the {ref['assembly']} "
                f"reference genome using BWA-MEM against the "
                f"FASTA specified in the run configuration. The Micro-C mapping "
                f"workflow used the BWA-MEM flags -5SP and -T0. "
                f"Aligned reads were parsed with pairtools using "
                f"a minimum mapping quality of "
                f"{processing['mapq']}, the "
                f"'{processing['walks_policy']}' walks policy and "
                f"a maximum inter-alignment gap of "
                f"{processing['max_inter_align_gap']} bp. "
                f"Read pairs were sorted and duplicate contact pairs "
                f"were removed with pairtools dedup. Downstream "
                f"analyses used the resulting non-duplicate pair "
                f"files."
            )
        )
    )

    sections.append(
        (
            "Contact matrices and balancing",
            (
                f"Contact matrices were generated with cooler. "
                f"The base contact-map resolution was "
                f"{fmt_bp(cmaps['base_resolution'])}, and "
                f"multi-resolution matrices were generated at "
                f"{join_bp(cmaps['resolutions'])}. Matrix balancing "
                f"excluded the first {cmaps['ignore_diags']} "
                f"diagonals. The primary chromosome view comprised "
                f"chr1-chr22 and chrX."
            )
        )
    )

    sections.append(
        (
            "Quality control and library complexity",
            (
                f"Quality-control criteria included a minimum of "
                f"{qc['deep_min_nodup_pairs']:,} non-duplicate pairs "
                f"for deep libraries and a minimum cis interaction "
                f"fraction at genomic separations >=1 kb of "
                f"{qc['min_cis_ge_1kb_fraction']:.2f}. "
                f"Genome-wide QC contact maps were evaluated at "
                f"{fmt_bp(qc['qc_map_resolution'])} resolution. "
                f"Library-complexity estimates were additionally "
                f"used to assess expected gains from deeper "
                f"sequencing."
            )
        )
    )

    sections.append(
        (
            "Reproducibility",
            (
                f"Contact-map reproducibility between samples was "
                f"assessed using HiCRep at "
                f"{fmt_bp(rep['resolution'])} resolution, with "
                f"smoothing parameter h={rep['hicrep_h']} and a "
                f"maximum genomic separation of "
                f"{fmt_bp(rep['hicrep_max_distance'])}. "
                f"Stratum-adjusted correlation coefficients were "
                f"calculated chromosome-wise and summarized across "
                f"chromosomes."
            )
        )
    )

    sections.append(
        (
            "A/B compartment analysis",
            (
                f"A/B compartment structure was evaluated using "
                f"cooltools eigs-cis. The primary compartment "
                f"analysis was performed at "
                f"{fmt_bp(compartments['primary_resolution'])} "
                f"resolution and eigenvectors were phased using "
                f"GC content. A configured absolute E1 threshold of "
                f"{compartments['primary_min_abs_e1']} was applied to screen for "
                f"strictly concordant sign switches (every control sample versus "
                f"every test sample). Such calls are descriptive, not statistically differential."
            )
        )
    )

    sections.append(
        (
            "Insulation and domain-boundary analysis",
            (
                f"Insulation scores were calculated with cooltools "
                f"at {fmt_bp(boundaries['resolution'])} resolution "
                f"using windows of "
                f"{join_bp(boundaries['windows'])}. "
                f"The primary insulation window was "
                f"{fmt_bp(boundaries['primary_window'])}, and "
                f"boundaries were identified using the "
                f"{boundaries['threshold_method']} thresholding "
                f"method. Boundary positions were clustered with a tolerance "
                f"of {fmt_bp(boundaries['cluster_tolerance_bp'])} and presence/absence "
                f"was compared across all samples of each group. These calls are "
                f"descriptive and are not FDR-significant differential-boundary tests."
            )
        )
    )

    sections.append(
        (
            "Genome-wide differential contact analysis",
            (
                f"Condition-blind cis contact features were "
                f"quantified at "
                f"{fmt_bp(contacts['resolution'])} resolution "
                f"between genomic separations of "
                f"{fmt_bp(contacts['min_distance_bp'])} and "
                f"{fmt_bp(contacts['max_distance_bp'])}. "
                f"Features were pre-filtered by requiring at least "
                f"{contacts['prefilter_min_total_count']} total "
                f"contacts and non-zero counts in at least "
                f"{contacts['prefilter_min_nonzero_samples']} "
                f"samples. Independent edgeR filtering used "
                f"min.count={contacts.get('edgeR_min_count', 5)} and "
                f"min.total.count={contacts.get('edgeR_min_total_count', 20)}. "
                f"Differential testing was performed with "
                f"the pipeline's edgeR implementation. Contacts "
                f"were classified as differential at "
                f"FDR < {contacts['differential_fdr']} and "
                f"|log2FC| >= "
                f"{contacts['differential_min_abs_log2fc']}."
            )
        )
    )

    sections.append(
        (
            "Loop detection",
            (
                f"Loop discovery was performed condition-blind from "
                f"a pooled contact map containing all {sample_count} samples. "
                f"Loops were called with Mustache (software version "
                f"recorded in the run provenance) at the "
                f"primary resolution of {fmt_bp(loop_res)}, using "
                f"a p-value threshold of "
                f"{loops['mustache_pvalue_threshold']} and a "
                f"sparsity threshold of "
                f"{loops['mustache_sparsity_threshold']}. "
                f"Only cis loop candidates between "
                f"{fmt_bp(loops['min_distance_bp'])} and "
                f"{fmt_bp(loops['max_distance_bp'])} were retained. "
                f"Candidate discovery was performed without using "
                f"condition labels."
            )
        )
    )

    sections.append(
        (
            "Differential loop analysis",
            (
                f"Each condition-blind loop candidate was quantified "
                f"separately in all {sample_count} samples using "
                f"raw contact counts from a {window_n}x{window_n}-bin "
                f"window centered on the loop pixel "
                f"({window_n * loop_res // 1000} kb by "
                f"{window_n * loop_res // 1000} kb at the primary "
                f"resolution). Balanced contact values were retained "
                f"for visualization and QC only. Differential loop "
                f"strength was tested using edgeR quasi-likelihood "
                f"negative-binomial models. Total full-depth contact "
                f"counts from each sample-specific cooler were used "
                f"as library-size exposure values and TMM factors "
                f"were used for between-sample composition "
                f"normalization. Low-count loops were filtered using "
                f"min.count={loops.get('edgeR_min_count', 5)} and "
                f"min.total.count="
                f"{loops.get('edgeR_min_total_count', 15)}. "
                f"Differential loops were defined as "
                f"FDR < {loops['differential_fdr']} and "
                f"|log2FC| >= "
                f"{loops['differential_min_abs_log2fc']}. "
                f"A stringent concordant core additionally required "
                f"every test sample to change in the same direction "
                f"as the group-level effect."
            )
        )
    )

    # The YAML may list optional resolutions that have not actually been run.
    # Record this *configuration* in provenance, but never describe a run as
    # completed in a manuscript without verified analysis outputs.
    secondary_resolutions = sorted({
        int(x) for x in sensitivity_res if int(x) != loop_res
    })
    if secondary_resolutions:
        sections.append(
            (
                "Configured optional loop resolutions (execution not verified)",
                (
                    f"The configuration lists {join_bp(secondary_resolutions)} "
                    f"as optional secondary loop resolutions. "
                    f"This Methods generator has not verified that those "
                    f"analyses were executed. They must not be reported as "
                    f"completed without separate run provenance."
                )
            )
        )

    if args.no_differential:
        sections = [(title, body) for title, body in sections
                    if title not in ("Genome-wide differential contact analysis", "Loop detection",
                                     "Differential loop analysis", "Configured optional loop resolutions (execution not verified)")]
        if sample_count == 1:
            sections = [(title, body) for title, body in sections
                        if title != "Reproducibility"]

    if cmaps.get("reuse_verified_existing_mcools", False):
        sections.append((
            "Legacy contact-map adoption (50-kb equivalence checked)",
            "Existing multi-resolution matrices were retained after a complete "
            "comparison of raw 50-kb pixels and bin coordinates against 10-kb "
            "coolers independently aggregated to 50 kb. Matching 50-kb raw "
            "counts do not by themselves certify other zoom levels or "
            "iterative-balancing weights. The supporting equivalence TSV is "
            "recorded under results/provenance."
        ))

    sections.append(
        (
            "Reproducibility and provenance",
            (
                f"The analysis was executed with the SnakeMicroC "
                f"workflow. The exact resolved configuration used "
                f"for this Methods text is archived with the results "
                f"(SHA256: {config_hash}). Pipeline git commit: "
                f"{commit}; repository status at report generation: "
                f"{git_status}. Software environments and version "
                f"information are retained in the analysis "
                f"provenance files."
            )
        )
    )

    # --------------------------------------------------------
    # Markdown
    # --------------------------------------------------------

    md_lines = [
        "# Micro-C analysis methods and technical provenance",
        "",
        (
            "This file was generated automatically from the "
            "resolved analysis configuration."
        ),
        "",
    ]

    for title, body in sections:
        md_lines.extend([
            f"## {title}",
            "",
            body,
            "",
        ])

    md_lines.extend([
        "## Parameter record",
        "",
        f"- Reference assembly: `{ref['assembly']}`",
        f"- MAPQ threshold: `{processing['mapq']}`",
        f"- Primary contact resolution: `{contacts['resolution']}` bp",
        f"- Primary loop resolution: `{loop_res}` bp",
        (
            "- Mustache p-value threshold: "
            f"`{loops['mustache_pvalue_threshold']}`"
        ),
        (
            "- Mustache sparsity threshold: "
            f"`{loops['mustache_sparsity_threshold']}`"
        ),
        (
            "- Loop differential FDR: "
            f"`{loops['differential_fdr']}`"
        ),
        (
            "- Loop minimum absolute log2FC: "
            f"`{loops['differential_min_abs_log2fc']}`"
        ),
        (
            "- Resolved configuration SHA256: "
            f"`{config_hash}`"
        ),
        "",
    ])

    md = "\n".join(md_lines)

    Path(args.output_md).parent.mkdir(
        parents=True,
        exist_ok=True
    )

    Path(args.output_md).write_text(md)

    # Paper-ready Methods version
    sec = dict(sections)

    mapping = sec["Read mapping and contact-pair processing"]
    mapping = re.sub(
        r"The Micro-C mapping workflow used.*?-T0\. ",
        "",
        mapping,
    )
    mapping = re.sub(
        r"Aligned reads were parsed with pairtools using.*?bp\. ",
        "Valid interaction pairs were identified with pairtools, retaining high-quality alignments with MAPQ >= {}. ".format(processing["mapq"]),
        mapping,
    )
    mapping = mapping.replace(
        "Read pairs were sorted and duplicate contact pairs were removed with pairtools dedup. ",
        "Duplicate contact pairs were removed before downstream analysis. ",
    )

    contact_maps = (
        "Contact matrices were generated with cooler and balanced using iterative "
        "correction. Multi-resolution contact maps were generated for downstream "
        "analyses, with resolution selected according to the genomic scale of the "
        "feature being evaluated."
    )

    reproducibility = (
        "Contact-map reproducibility between samples was assessed on autosomes "
        "using HiCRep at {} resolution, considering interactions up to {} genomic "
        "distance.".format(
            fmt_bp(rep["resolution"]),
            fmt_bp(rep["hicrep_max_distance"]),
        )
    )

    compartment = (
        "A/B compartment organization was evaluated at {} resolution using "
        "eigenvector decomposition implemented in cooltools. The first "
        "eigenvector was oriented using GC content to assign chromatin "
        "compartment state.".format(
            fmt_bp(compartments["primary_resolution"])
        )
    )

    insulation = (
        "Local chromatin-domain organization was assessed using insulation-score "
        "analysis at {} resolution with {} windows. The {} window was used as "
        "the primary scale for boundary analysis, and domain boundaries were "
        "identified using the {} thresholding approach.".format(
            fmt_bp(boundaries["resolution"]),
            join_bp(boundaries["windows"]),
            fmt_bp(boundaries["primary_window"]),
            boundaries["threshold_method"],
        )
    )

    differential_contacts = (
        "Genome-wide differential chromatin interactions were evaluated on "
        "autosomes at {} resolution. Cis interactions spanning {} to {} were "
        "quantified independently in each sample, and low-count interactions "
        "were removed before statistical testing. Differential contacts between "
        "test and control samples were identified using edgeR. Contacts "
        "were considered differential at FDR < {} and |log2FC| >= {}.".format(
            fmt_bp(contacts["resolution"]),
            fmt_bp(contacts["min_distance_bp"]),
            fmt_bp(contacts["max_distance_bp"]),
            contacts["differential_fdr"],
            contacts["differential_min_abs_log2fc"],
        )
    )

    loop_detection = sec.get("Loop detection", "")

    differential_loops = sec.get("Differential loop analysis", "")
    differential_loops = re.sub(
        r"Total full-depth contact counts.*?normalization\. ",
        "TMM normalization was used to account for between-sample composition differences. ",
        differential_loops,
    )
    differential_loops = re.sub(
        r"Low-count loops were filtered.*?\. ",
        "Low-count loops were removed before statistical testing. ",
        differential_loops,
    )

    paper_sections = [
        (
            "Micro-C data processing and contact map generation",
            sec["Experimental design and sequencing"] + " " + mapping + " " + contact_maps,
        ),
        (
            "Reproducibility, compartment and insulation analyses",
            (reproducibility + " " if sample_count > 1 else "") + compartment + " " + insulation,
        ),
    ]
    if not args.no_differential:
        paper_sections.extend([
            ("Differential chromatin contacts", differential_contacts),
            ("Chromatin loop detection and differential loop analysis", loop_detection + " " + differential_loops),
        ])

    # Do not add a sensitivity-analysis paragraph to manuscript Methods:
    # a configured resolution is not proof of a completed analysis.

    txt_lines = []

    for title, body in paper_sections:
        txt_lines.extend([
            title,
            "",
            body,
            "",
        ])

    Path(args.output_txt).write_text(
        "\n".join(txt_lines)
    )

    print(
        f"Methods Markdown: {args.output_md}"
    )

    print(
        f"Methods text: {args.output_txt}"
    )


if __name__ == "__main__":
    main()
