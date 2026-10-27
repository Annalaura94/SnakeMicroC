"""Regression checks for pairtools scaling's named-column --view format."""
from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from prepare_pairtools_scaling_view import HEADER, make_view


class PairtoolsScalingViewTests(unittest.TestCase):
    def test_headerless_bed4_converts_and_preserves_input(self):
        with tempfile.TemporaryDirectory() as folder:
            base = Path(folder)
            inp = base / "main.view.tsv"
            out = base / "pairtools.view.tsv"
            original = "chr1\t0\t10000\tchr1\nchr2\t0\t20000\tchr2\n"
            inp.write_text(original)
            self.assertEqual(make_view(inp, out), 2)
            self.assertEqual(inp.read_text(), original)
            with out.open(newline="") as handle:
                rows = list(csv.reader(handle, delimiter="\t"))
            self.assertEqual(rows[0], list(HEADER))
            self.assertEqual(rows[1:], [
                ["chr1", "0", "10000", "chr1"],
                ["chr2", "0", "20000", "chr2"],
            ])

    def test_invalid_missing_name_and_duplicates_are_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            base = Path(folder)
            inp = base / "main.tsv"
            out = base / "pairtools.tsv"
            for content in (
                "chr1\t0\t5000\n",
                "chr1\t0\t5000\tchr1\nchr2\t0\t8000\tchr1\n",
                "chr1\t5000\t1000\tchr1\n",
            ):
                inp.write_text(content)
                with self.assertRaises(ValueError):
                    make_view(inp, out)
                self.assertFalse(out.exists())

    def test_snakefile_uses_independent_scaling_view_and_isolated_python(self):
        rules = (ROOT / "rules/Snakefile").read_text()
        self.assertIn('SCALING_VIEW = f"{RESULTS}/reference/main_chromosomes.pairtools.view.tsv"', rules)
        self.assertIn("rule prepare_pairtools_scaling_view:", rules)
        self.assertIn("view=SCALING_VIEW", rules)
        self.assertIn("export PYTHONNOUSERSITE=1", rules)
        self.assertIn("main_view=MAIN_VIEW", rules)
        self.assertIn("--view {input.view}", rules)


if __name__ == "__main__":
    unittest.main()
