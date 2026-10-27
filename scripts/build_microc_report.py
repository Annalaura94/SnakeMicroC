#!/usr/bin/env python3
"""Generate one navigable technical report from completed SnakeMicroC outputs.

The report is strictly evidence-based: data tables must exist for a section to
be marked completed. It does not infer biological significance from file names.
The input paths are relative to the Snakemake experiment workspace.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import html
import os
from pathlib import Path
from urllib.parse import quote

import yaml


def load_table(path):
    with Path(path).open(encoding="utf8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def metrics(path):
    return {str(x.get("metric", "")): str(x.get("value", "")) for x in load_table(path)}


def structural_direction_summary(path):
    """Descriptive sample-direction concordance, never a significance test."""
    opener = gzip.open if str(path).endswith(".gz") else open
    counts = {"tested_bins": 0, "concordant_bins": 0,
              "concordant_up": 0, "concordant_down": 0}
    with opener(path, "rt", encoding="utf8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            counts["tested_bins"] += 1
            direction = str(row.get("consensus_direction", "")).lower()
            if direction == "up":
                counts["concordant_up"] += 1
                counts["concordant_bins"] += 1
            elif direction == "down":
                counts["concordant_down"] += 1
                counts["concordant_bins"] += 1
            elif str(row.get("all_test_same_direction", "")).lower() == "true":
                counts["concordant_bins"] += 1
    return counts


def html_table(rows, columns=None, limit=200):
    if not rows:
        return "<p class='muted'>No rows available.</p>"
    columns = columns or list(rows[0])
    esc = html.escape
    header = "".join(f"<th>{esc(str(col))}</th>" for col in columns)
    body = "\n".join("<tr>" + "".join(f"<td>{esc(str(row.get(c, '')))}</td>" for c in columns) + "</tr>" for row in rows[:limit])
    tail = f"<p class='muted'>Showing the first {limit:,} rows.</p>" if len(rows) > limit else ""
    return f"<div class='table-scroll'><table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table></div>{tail}"


def read_text(path):
    p = Path(path)
    return p.read_text(encoding="utf8", errors="replace") if p.is_file() else ""


def file_url(file_path, report_dir):
    rel = os.path.relpath(file_path, report_dir)
    return quote(rel.replace(os.sep, "/"), safe="/-._")


def render_report(config, results, report_dir, differential):
    samples = config["samples"]
    condition_counts = {}
    for entry in samples.values():
        cond = entry["condition"]
        condition_counts[cond] = condition_counts.get(cond, 0) + 1
    comp = config.get("comparison") or {}
    ctrl, test = comp.get("control_condition", ""), comp.get("test_condition", "")
    sample_label = ", ".join(f"{s} ({v['condition']})" for s, v in samples.items())
    esc = html.escape
    sections = []
    index = []
    summary = []
    expected = [
        ("QC", "qc/qc_summary.tsv", "QC metrics per sample", True),
        ("QC", "qc/qc_status.txt", "QC pass / fail checks", True),
        ("QC", "qc/library_complexity.tsv", "Complexity and sequencing depth", True),
        ("QC", "qc/scaling_curves.tsv", "Contact decay curve", True),
        ("QC", "qc/contact_map_correlation_1Mb.tsv", "Inter-sample correlation", True),
        ("Matrices", "mcool/{sample}.mcool", "Multi-resolution contact matrix", True),
        ("Structure", "reproducibility/hicrep_SCC_matrix.tsv", "HiCRep SCC", len(samples) > 1),
        ("Structure", "structural/compartments_clone_effects.tsv.gz", "Compartment comparison", differential),
        ("Structure", "structural/insulation_clone_effects.tsv.gz", "Insulation comparison", differential),
        ("Contacts", "contacts/differential_contacts50_summary.tsv", "Differential contact summary", differential),
        ("Loops", f"loops_mustache/{config['loops']['primary_resolution']}/differential_loop_summary.tsv", "Differential loop summary", differential),
        ("Methods", "report/MicroC_methods_paper.txt", "Paper-ready methods", True),
        ("Provenance", "provenance/software_versions.tsv", "Resolved software versions", True),
        ("Provenance", "provenance/resolved_config.yaml", "Resolved configuration", True),
        ("Provenance", "provenance/mcool_50kb_pixel_equivalence.tsv", "Full 50-kb raw-pixel equivalence for adopted legacy matrices",
         bool(config.get("contact_maps", {}).get("reuse_verified_existing_mcools", False))),
    ]
    # Convenient links to canonical machine-readable results. These optional
    # artifacts do not block report creation if the relevant module is disabled.
    expected.extend([
        ("QC", "qc/scaling_curves.pdf", "Distance-decay plot", True),
        ("QC", "qc/library_complexity.pdf", "Complexity figure", True),
        ("Structure", "structural/compartments/{sample}." + str(config['compartments']['primary_resolution']) + ".cis.vecs.tsv", "Per-sample E1", True),
        ("Structure", "structural/insulation/{sample}." + str(config['boundaries']['resolution']) + "." + str(config['boundaries']['primary_window']) + ".tsv", "Per-sample insulation", True),
        ("Contacts", "contacts/differential_contacts50_all.tsv.gz", "All tested differential contacts", differential),
        ("Contacts", "contacts/differential_contacts50_core.tsv", "High-confidence contact candidates", differential),
        ("Contacts", "contacts/contact50_prefilter_summary.tsv", "Distance and count prefilter summary", differential),
        ("Contacts", "contacts/EXP_independent_contacts50_volcano.pdf", "Contact volcano plot", differential),
        ("Loops", f"loops_mustache/{config['loops']['primary_resolution']}/differential_loops_all.tsv.gz", "All tested loops", differential),
        ("Loops", f"loops_mustache/{config['loops']['primary_resolution']}/differential_loops_core.tsv", "Concordant loop candidates", differential),
        ("Loops", f"loops_mustache/{config['loops']['primary_resolution']}/MicroC_loop_analysis_report.html", "Full loop report", differential),
        ("Loops", f"loops_mustache/{config['loops']['primary_resolution']}/pooled_loop_candidates_counts.tsv.gz", "Raw window counts", differential),
    ])
    for group, pattern, desc, required in expected:
        paths = [pattern.format(sample=s) for s in samples] if "{sample}" in pattern else [pattern]
        for rp in paths:
            filepath = results / rp
            state = "ready" if filepath.is_file() and filepath.stat().st_size > 0 else ("missing" if required else "not-applicable")
            index.append({"module": group, "description": desc, "status": state, "path": str(filepath) if state == "ready" else rp})
    missing = [x for x in index if x["status"] == "missing"]
    if missing:
        raise RuntimeError("Required report inputs missing: " + "; ".join(x["path"] for x in missing))

    def artifact(rel, label, preview=False):
        path = results / rel
        if not path.is_file():
            return ""
        href = file_url(path, report_dir)
        img = f"<object data='{esc(href)}' type='application/pdf' class='pdf-preview'><p>Open PDF using the link above.</p></object>" if preview and path.suffix == '.pdf' else ""
        return f"<div class='artifact'><a href='{esc(href)}'>{esc(label)}</a>{img}</div>"

    def metric_summary(rel, section):
        rows = metrics(results / rel)
        for name, value in rows.items():
            summary.append({"section": section, "metric": name, "value": value})
        return html_table([{"metric": k, "value": v} for k, v in rows.items()])

    qc = load_table(results / "qc/qc_summary.tsv")
    qc_cols = [x for x in ["sample", "condition", "total_read_pairs", "no_dup_read_pairs", "mapped_pct_total", "duplicates_pct_total", "cis_ge_1kb_pct_nodup", "dovetail_deep_qc_pass"] if x in qc[0]] if qc else []
    for row in qc:
        summary.append({"section": "QC", "metric": f"{row['sample']}.no_dup_read_pairs", "value": str(row.get('no_dup_read_pairs', ''))})
        summary.append({"section": "QC", "metric": f"{row['sample']}.deep_qc_pass", "value": str(row.get('dovetail_deep_qc_pass', ''))})
    sections.append("<section id='quality'><h2>1. Read-pair quality control</h2>"
                    "<p>These metrics are extracted from pairtools statistics. A failed QC threshold is shown as measured, not masked.</p>"
                    + html_table(qc, qc_cols)
                    + artifact("qc/qc_summary.pdf", "Full QC plots (PDF)", preview=True)
                    + artifact("qc/library_complexity.pdf", "Library-complexity estimates (PDF)", preview=False)
                    + artifact("qc/scaling_curves.pdf", "Contact decay P(s) (PDF)", preview=False)
                    + artifact("qc/contact_maps_1Mb.pdf", "One-megabase contact maps (PDF)", preview=False)
                    + artifact("qc/contact_map_correlation_1Mb.pdf", "Inter-sample contact map correlation (PDF)", preview=False)
                    + "</section>")
    sections.append("<section id='structure'><h2>2. Reproducibility and 3D structure</h2>"
                    "<p>Contact matrices, HiCRep, A/B compartments and insulation were evaluated as separate structural layers.</p>")
    rep = results / "reproducibility/hicrep_SCC_matrix.tsv"
    if rep.is_file() and len(samples) > 1:
        sections.append("<h3>HiCRep SCC</h3>" + html_table(load_table(rep)))
    else:
        sections.append("<p>HiCRep not applicable for a single sample.</p>")
    if differential:
        sections.append("<h3>Descriptive compartment and insulation effects</h3>"
                        "<p>Direction concordance of E1 and insulation-score differences across test samples. These counts are not statistical discoveries, compartment switches or boundary-gain calls.</p>")
        for label, rel in (("A/B E1", "structural/compartments_clone_effects.tsv.gz"),
                           ("Insulation", "structural/insulation_clone_effects.tsv.gz")):
            counts = structural_direction_summary(results / rel)
            for key, value in counts.items():
                summary.append({"section": label, "metric": key, "value": value})
            sections.append("<h3>" + esc(label) + "</h3>"
                            + html_table([{"metric": k, "value": v} for k, v in counts.items()]))
    sections.append("</section>")
    if differential:
        sections.append("<section id='differential'><h2>3. Differential interactions</h2>"
                        "<p><strong>Statistical caution:</strong> sample count and independent biological replication limit power. Differential calls depend on filtering and normalization."
                        " Chromatin loop contact counts must not be treated as direct proof of focal loop formation or loss.</p>"
                        "<h3>Contacts at configured resolution</h3>"
                        + metric_summary("contacts/differential_contacts50_summary.tsv", "Differential contacts")
                        + artifact("contacts/EXP_independent_contacts50_volcano.pdf", "Differential contacts volcano (PDF)", preview=True)
                        + artifact("contacts/contact50_prefilter_summary.tsv", "Condition-blind contact prefilter table")
                        + "<h3>Mustache pooled loops</h3>"
                        + metric_summary(f"loops_mustache/{config['loops']['primary_resolution']}/differential_loop_summary.tsv", "Differential loops")
                        + artifact(f"loops_mustache/{config['loops']['primary_resolution']}/EXP_differential_loops_volcano.pdf", "Differential loops volcano (PDF)")
                        + artifact(f"loops_mustache/{config['loops']['primary_resolution']}/EXP_differential_loops_QC.pdf", "Differential loops QC (PDF)", preview=True)
                        + artifact(f"loops_mustache/{config['loops']['primary_resolution']}/MicroC_loop_analysis_report.html", "Detailed loop report (HTML)")
                        + "</section>")
    else:
        sections.append("<section id='differential'><h2>3. Differential interactions</h2><p><strong>Not run:</strong> fewer than two samples per group, or differential analysis was disabled in the configuration. This run still provides QC and per-sample structural results.</p></section>")
    methods = read_text(results / "report/MicroC_methods_paper.txt")
    sections.append("<section id='methods'><h2>4. Methods</h2><p>Automatically generated from the resolved experiment configuration and the completed workflow targets. Verify wording before manuscript submission.</p>"
                    "<pre class='methods'>" + esc(methods) + "</pre>"
                    + artifact("report/MicroC_methods_provenance.md", "Methods with technical parameters (Markdown)")
                    + "</section>")
    def index_table():
        header = "<thead><tr><th>Module</th><th>Purpose</th><th>Status</th><th>File</th></tr></thead>"
        rows = []
        for item in index:
            path = Path(item["path"])
            if item["status"] == "ready":
                cell = f"<a href='{esc(file_url(path, report_dir))}'>{esc(str(path.relative_to(results)))}</a>"
            else:
                cell = esc(item["path"])
            rows.append("<tr>" + "".join(f"<td>{esc(str(item[key]))}</td>" for key in ("module", "description", "status")) + f"<td>{cell}</td></tr>")
        return "<div class='table-scroll'><table>" + header + "<tbody>" + "".join(rows) + "</tbody></table></div>"

    if config.get("contact_maps", {}).get("reuse_verified_existing_mcools", False):
        summary.append({"section": "Provenance", "metric": "mcool_matrix_origin", "value": "audited-existing"})
        adoption_note = ("<p><strong>Audited legacy matrices:</strong> existing .mcool files were retained "
                         "following a complete raw-pixel comparison at 50 kb against independently "
                         "aggregated 10-kb coolers. Other resolutions and balancing weights were "
                         "not independently certified by that comparison.</p>")
    else:
        summary.append({"section": "Provenance", "metric": "mcool_matrix_origin", "value": "workflow-generated"})
        adoption_note = ""
    sections.append("<section id='files'><h2>5. Canonical outputs and software</h2>"
                    "<p>Primary scientific results remain in structured TSV/TSV.GZ and cooler files. This index points to canonical outputs, not exploratory backups.</p>"
                    + adoption_note + index_table()
                    + artifact("provenance/software_versions.tsv", "Software versions (TSV)")
                    + artifact("provenance/resolved_config.yaml", "Resolved configuration (YAML)")
                    + "</section>")
    report_title = "SnakeMicroC | experiment report"
    style = """
    :root{font-family:Inter,system-ui,-apple-system,sans-serif;color:#162337;background:#f5f7fa}
    body{max-width:1120px;margin:0 auto;padding:36px 24px;line-height:1.55}
    h1{font-size:2.3rem;letter-spacing:-.04em;margin-bottom:6px}
    h2{font-size:1.5rem;margin-top:0}h3{font-size:1.15rem}
    .muted{color:#59687a}.lead{font-size:1.05rem;color:#415166}
    .card,section{background:white;padding:26px;border:1px solid #dce3ec;border-radius:12px;margin:22px 0}
    .card{background:#ebf2fa}.table-scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:.9rem}
    th,td{padding:10px;text-align:left;border-bottom:1px solid #e4e9ee;white-space:nowrap}
    th{background:#f3f6fa;font-size:.8rem;text-transform:uppercase;letter-spacing:.04em;color:#43566e}
    a{color:#135ea8;text-decoration:none;font-weight:600}a:hover{text-decoration:underline}
    .artifact{margin:12px 0}.pdf-preview{display:block;width:100%;min-height:420px;margin-top:10px;border:1px solid #e2e8f0}
    pre.methods{white-space:pre-wrap;font-size:.9rem;background:#f7f9fc;padding:18px;line-height:1.65}
    nav a{display:inline-block;margin-right:20px;padding:8px 0}
    footer{font-size:.84rem;color:#63758a;margin-top:42px}
    """
    index_html = "<nav>" + "".join(f"<a href='#{key}'>{title}</a>" for key, title in [
        ("quality", "QC"), ("structure", "Structure"), ("differential", "Differential"),
        ("methods", "Methods"), ("files", "Files")]) + "</nav>"
    heading = (f"<header><h1>{esc(report_title)}</h1><p class='lead'>{len(samples)} samples · {esc(sample_label)}</p>"
               f"<div class='card'><strong>Comparison:</strong> {esc((test + ' versus ' + ctrl) if differential else 'not applicable')}"
               f"<p>Reference: {esc(str(config['reference']['assembly']))} · primary loops: {esc(str(config['loops']['primary_resolution']))} bp"
               f" · contacts: {esc(str(config['contacts']['resolution']))} bp</p></div></header>")
    html_doc = ("<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width'>"
                f"<title>{esc(report_title)}</title><style>{style}</style></head><body>"
                + heading + index_html + "\n".join(sections)
                + "<footer>Generated by SnakeMicroC. This report summarizes pipeline measurements; biological interpretation and external validation remain the investigator's responsibility.</footer></body></html>")
    return html_doc, summary, index, condition_counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--results", required=True)
    parser.add_argument("--report-dir", required=True)
    parser.add_argument("--differential", choices=("yes", "no"), required=True)
    args = parser.parse_args()
    cfg = yaml.safe_load(Path(args.config).read_text(encoding="utf8"))
    results = Path(args.results).resolve()
    report_dir = Path(args.report_dir).resolve()
    report_dir.mkdir(parents=True, exist_ok=True)
    html_text, summary, index, conditions = render_report(cfg, results, report_dir, args.differential == "yes")
    (report_dir / "MicroC_report.html").write_text(html_text, encoding="utf8")
    for name, rows, cols in [
        ("MicroC_summary.tsv", summary, ("section", "metric", "value")),
        ("MicroC_output_index.tsv", index, ("module", "description", "status", "path")),
    ]:
        with (report_dir / name).open("w", encoding="utf8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    print(f"Generated: {report_dir / 'MicroC_report.html'}; groups: {conditions}")


if __name__ == "__main__":
    main()
