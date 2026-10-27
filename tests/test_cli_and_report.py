"""Fast regression tests for CLI, configuration and integrated reporting.

These tests use synthetic TSVs and fake Snakemake, never Micro-C raw data.
"""
from __future__ import annotations

import csv
import gzip
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from validate_config import ConfigurationError, load_and_validate


def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def new_config(tmp, n_control=2, n_test=2):
    cfg = yaml.safe_load((ROOT / "config/config.template.yaml").read_text())
    cfg["tmpdir"] = str(tmp / "scratch")
    cfg["reference"]["fasta"] = str(tmp / "reference.fa")
    cfg["reference"]["fai"] = str(tmp / "reference.fa.fai")
    (tmp / "reference.fa").write_text(">chr1\nACGT\n")
    (tmp / "reference.fa.fai").write_text("chr1\t4\t6\t4\t5\n")
    samples = {}
    for cond, n in (("CONTROL", n_control), ("TEST", n_test)):
        for i in range(1, n+1):
            name = f"{cond}_{i}"
            r1, r2 = tmp / f"{name}_R1.fastq.gz", tmp / f"{name}_R2.fastq.gz"
            for path in (r1, r2):
                with gzip.open(path, "wt") as fh:
                    fh.write("@read\nACGT\n+\nIIII\n")
            samples[name] = {"condition": cond, "clone": name, "r1": str(r1), "r2": str(r2)}
    cfg["samples"] = samples
    cfg["contacts"]["prefilter_min_nonzero_samples"] = min(len(samples), 2)
    cfgpath = tmp / "config.yaml"
    cfgpath.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return cfg, cfgpath


class ConfigTests(unittest.TestCase):
    def test_new_dataset_preflight(self):
        with tempfile.TemporaryDirectory() as name:
            tmp = Path(name)
            _, config = new_config(tmp, 3, 2)
            info = load_and_validate(config, check_paths=True)
            self.assertEqual(info["samples"], 5)
            self.assertTrue(info["differential_enabled"])

    def test_single_sample_skips_differential(self):
        with tempfile.TemporaryDirectory() as name:
            tmp = Path(name)
            cfg, path = new_config(tmp, 1, 0)
            cfg.pop("comparison")
            cfg["analysis"]["require_differential"] = False
            path.write_text(yaml.safe_dump(cfg))
            info = load_and_validate(path, check_paths=True)
            self.assertFalse(info["differential_enabled"])
            self.assertEqual(info["samples"], 1)

    def test_bad_sample_id_and_wrong_resolution_rejected(self):
        with tempfile.TemporaryDirectory() as name:
            tmp = Path(name)
            cfg, path = new_config(tmp, 2, 2)
            cfg["samples"]["Bad sample"] = cfg["samples"].pop("CONTROL_1")
            path.write_text(yaml.safe_dump(cfg))
            with self.assertRaises(ConfigurationError):
                load_and_validate(path)
            cfg, path = new_config(tmp, 2, 2)
            cfg["contact_maps"]["resolutions"].remove(100000)
            path.write_text(yaml.safe_dump(cfg))
            with self.assertRaisesRegex(ConfigurationError, "resolutions"):
                load_and_validate(path)

    def test_cli_passes_correct_snakemake_target(self):
        with tempfile.TemporaryDirectory() as name:
            tmp = Path(name)
            _, cfg = new_config(tmp)
            fake = tmp / "fake-snakemake"
            log = tmp / "args.json"
            fake.write_text("#!/usr/bin/env python3\nimport json,os,sys\nopen(os.environ['CALL_LOG'],'w').write(json.dumps(sys.argv[1:]))\n")
            fake.chmod(0o755)
            for subcommand, target in (("qc", "qc"), ("run", "all"), ("report", "report")):
                result = subprocess.run([
                    sys.executable, str(ROOT / "bin/microc"), subcommand,
                    "--config", str(cfg), "--snakemake", str(fake), "--dry-run", "--cores", "3",
                ], env={**os.environ, "CALL_LOG": str(log)}, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                args = json.loads(log.read_text())
                self.assertEqual(args[-1], target)
                self.assertIn("--use-conda", args)
                self.assertIn("--dry-run", args)
                self.assertIn("--directory", args)


class ReportTests(unittest.TestCase):
    def make_results(self, tmp, differential):
        cfg, path = new_config(tmp, 2 if differential else 1, 2 if differential else 0)
        if not differential:
            cfg.pop("comparison")
            path.write_text(yaml.safe_dump(cfg))
        results = tmp / "results"
        write_tsv(results / "qc/qc_summary.tsv", [
            {"sample": k, "condition": v["condition"], "no_dup_read_pairs": 100000, "dovetail_deep_qc_pass": "True"}
            for k, v in cfg["samples"].items()
        ])
        (results / "qc/qc_status.txt").write_text("CONTROL_1\tPASS\n")
        for rel in ("qc/scaling_curves.pdf", "qc/library_complexity.pdf"):
            p = results / rel
            p.write_bytes(b"fake pdf fixture")
        for sample in cfg["samples"]:
            comp = results / f"structural/compartments/{sample}.100000.cis.vecs.tsv"
            ins = results / f"structural/insulation/{sample}.10000.500000.tsv"
            for p in (comp, ins):
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text("chrom\tstart\tend\nchr1\t100\t200\n")
        write_tsv(results / "qc/library_complexity.tsv", [{"sample":"CONTROL_1", "depth": 100}])
        write_tsv(results / "qc/scaling_curves.tsv", [{"sample":"CONTROL_1", "distance_bp": 1000}])
        write_tsv(results / "qc/contact_map_correlation_1Mb.tsv", [{"sample":"CONTROL_1", "CONTROL_1": 1}])
        for sample in cfg["samples"]:
            p = results / f"mcool/{sample}.mcool"
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(b"synthetic not actual cool")
        if len(cfg["samples"]) > 1:
            write_tsv(results / "reproducibility/hicrep_SCC_matrix.tsv", [{"sample":"CONTROL_1", "CONTROL_1":1}])
        if differential:
            for rel in ("structural/compartments_clone_effects.tsv.gz", "structural/insulation_clone_effects.tsv.gz"):
                p=results/rel; p.parent.mkdir(parents=True,exist_ok=True)
                with gzip.open(p, "wt") as fh:
                    fh.write("chrom\tstart\tend\tconsensus_direction\nchr1\t100\t200\tup\n")
            for rel in (
                "contacts/differential_contacts50_all.tsv.gz",
                "contacts/differential_contacts50_core.tsv",
                "contacts/contact50_prefilter_summary.tsv",
                "contacts/EXP_independent_contacts50_volcano.pdf",
                "loops_mustache/10000/differential_loops_all.tsv.gz",
                "loops_mustache/10000/differential_loops_core.tsv",
                "loops_mustache/10000/MicroC_loop_analysis_report.html",
                "loops_mustache/10000/pooled_loop_candidates_counts.tsv.gz",
            ):
                p = results / rel
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_bytes(b"synthetic")
            write_tsv(results / "contacts/differential_contacts50_summary.tsv", [{"metric":"tested_after_filter", "value":122}, {"metric":"FDR_significant","value":4}])
            write_tsv(results / "loops_mustache/10000/differential_loop_summary.tsv", [{"metric":"condition_blind_candidates", "value":9535}, {"metric":"FDR_significant","value":1}])
        (results / "report").mkdir(parents=True)
        (results / "report/MicroC_methods_paper.txt").write_text("This is a synthetic METHODS section.")
        (results / "provenance").mkdir(parents=True)
        write_tsv(results / "provenance/software_versions.tsv", [{"tool":"python", "version":"3"}])
        (results / "provenance/resolved_config.yaml").write_text(path.read_text())
        return path, results

    def test_report_differential_and_single_sample(self):
        for differential in (False, True):
            with self.subTest(differential=differential), tempfile.TemporaryDirectory() as name:
                tmp=Path(name)
                cfg, results=self.make_results(tmp, differential)
                result=subprocess.run([
                    sys.executable, str(ROOT/"scripts/build_microc_report.py"),
                    "--config",str(cfg),"--results",str(results),"--report-dir",str(results/"report"),
                    "--differential","yes" if differential else "no"
                ],capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stderr)
                report=(results/"report/MicroC_report.html").read_text()
                self.assertIn("This is a synthetic METHODS",report)
                self.assertIn("CONTRO", report)
                if differential:
                    self.assertIn("FDR_significant", report)
                if not differential:
                    self.assertIn("Not run:",report)
                    self.assertNotIn("Mustache pooled loops",report)
                else:
                    self.assertIn("Mustache pooled loops",report)
                self.assertTrue((results/"report/MicroC_output_index.tsv").is_file())

    def test_methods_single_sample_no_claims(self):
        with tempfile.TemporaryDirectory() as name:
            tmp=Path(name)
            cfg,path=new_config(tmp,1,0)
            cfg.pop("comparison")
            path.write_text(yaml.safe_dump(cfg))
            paper=tmp/"methods.txt"
            md=tmp/"methods.md"
            out=subprocess.run([
                sys.executable,str(ROOT/"scripts/build_microc_methods.py"),
                "--config",str(path),"--pipeline-dir",str(ROOT),
                "--no-differential", "--output-md",str(md),"--output-txt",str(paper),
            ],capture_output=True,text=True)
            self.assertEqual(out.returncode,0,out.stderr)
            self.assertNotIn("Differential chromatin contacts",paper.read_text())
            self.assertNotIn("Loop detection",md.read_text())
            self.assertIn("1 libraries",md.read_text())


if __name__ == "__main__":
    unittest.main()
