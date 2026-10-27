"""Regression checks for strict, optional legacy .mcool adoption.

These tests never touch real data and do not depend on cooler/Snakemake.
"""
from __future__ import annotations

import csv
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
from verified_mcool import AdoptionError, verify_existing_mcools


class MatrixAdoptionTests(unittest.TestCase):
    def fixture(self, tmp):
        results = tmp / "results"
        samples = ["S1", "S2"]
        cfg = yaml.safe_load((ROOT / "config/config.template.yaml").read_text())
        cfg["analysis"]["require_differential"] = False
        cfg["samples"] = {s: {"condition": "C", "clone": s, "r1": f"/fake/{s}_R1.fq.gz", "r2": f"/fake/{s}_R2.fq.gz"} for s in samples}
        manifest = results / "provenance/mcool_50kb_pixel_equivalence.tsv"
        manifest.parent.mkdir(parents=True, exist_ok=True)
        for s in samples:
            for rel in (f"mcool/{s}.mcool", f"loops_mustache/10000/cool/{s}.cool"):
                x = results / rel
                x.parent.mkdir(parents=True, exist_ok=True)
                x.write_bytes(b"fake-matrix")
        with manifest.open("w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(("sample", "bins_match", "pixels_match", "status"))
            for s in samples:
                writer.writerow((s, "YES", "YES", "PASS"))
        cfg["contact_maps"]["reuse_verified_existing_mcools"] = True
        cfg["contact_maps"]["existing_mcool_equivalence_tsv"] = str(manifest)
        return cfg, manifest

    def test_valid_audit_adopts_without_reprocessing(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            cfg, manifest = self.fixture(tmp)
            self.assertEqual(verify_existing_mcools(cfg, tmp), manifest)
            cfg["contact_maps"]["reuse_verified_existing_mcools"] = False
            self.assertIsNone(verify_existing_mcools(cfg, tmp))

    def test_rejects_failed_missing_or_stale_audit(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            cfg, manifest = self.fixture(tmp)
            text = manifest.read_text().replace("S2\tYES\tYES\tPASS", "S2\tYES\tNO\tFAIL")
            manifest.write_text(text)
            with self.assertRaisesRegex(AdoptionError, "Every configured"):
                verify_existing_mcools(cfg, tmp)
            self.fixture(tmp)
            extra = tmp / "results/mcool/S1.mcool"
            os.utime(extra, (manifest.stat().st_mtime + 20, manifest.stat().st_mtime + 20))
            with self.assertRaisesRegex(AdoptionError, "newer than equivalence"):
                verify_existing_mcools(cfg, tmp)
            extra.unlink()
            with self.assertRaisesRegex(AdoptionError, "Missing/nonempty matrix"):
                verify_existing_mcools(cfg, tmp)

    def test_cli_rerun_trigger_modes(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            cfg, manifest = self.fixture(tmp)
            cfg["contact_maps"]["reuse_verified_existing_mcools"] = False
            cfg_file = tmp / "cfg.yaml"
            cfg_file.write_text(yaml.safe_dump(cfg))
            fake = tmp / "fake-snakemake"
            log = tmp / "args.json"
            fake.write_text("#!/usr/bin/env python3\nimport json,os,sys\nopen(os.environ['CALL_LOG'],'w').write(json.dumps(sys.argv[1:]))\n")
            fake.chmod(0o755)
            for mode in ("full", "mtime", "code-and-params"):
                args = [sys.executable, str(ROOT / "bin/microc"), "run", "--config", str(cfg_file),
                        "--workdir", str(tmp), "--skip-path-checks", "--snakemake", str(fake),
                        "--rerun-triggers", mode, "--dry-run"]
                result = subprocess.run(args, capture_output=True, text=True, env={**os.environ, "CALL_LOG": str(log)})
                self.assertEqual(result.returncode, 0, result.stderr)
                observed = json.loads(log.read_text())
                at = observed.index("--rerun-triggers")
                selected = observed[at+1:]
                for opt in ("--conda-frontend", "--conda-prefix", "--dry-run", "--printshellcmds", "all"):
                    if opt in selected:
                        selected = selected[:selected.index(opt)]
                self.assertEqual(set(selected), {
                    "full": {"mtime", "input", "params", "code", "software-env"},
                    "mtime": {"mtime"},
                    "code-and-params": {"input", "params", "code", "software-env"},
                }[mode])

    def test_rules_share_existing_10kb_cooler_and_mcool_for_qc(self):
        source = (ROOT / "rules/Snakefile").read_text()
        self.assertNotIn("rule make_base_cool:", source)
        self.assertNotIn("rule make_qc_cool:", source)
        self.assertIn("if not REUSE_VERIFIED_MCOOLS:", source)
        self.assertIn("if INS_RES != LOOP_RES:", source)
        self.assertIn('cool=INS_COOL', source)
        self.assertIn('cool_uris=QC_MCOOL_URIS', source)
        self.assertIn('configuration_sha256=CONFIG_FINGERPRINT', source)


if __name__ == "__main__":
    unittest.main()
