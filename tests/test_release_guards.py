"""Small regression tests for the first SnakeMicroC consolidation patch.

Run from the repository root:
    python -m unittest discover -s tests -p 'test_release_guards.py' -v

No FASTQ or cluster job execution is necessary. Requires the YAML, NumPy
and pandas dependencies already used by SnakeMicroC.
"""

from __future__ import annotations

import gzip
import importlib.util
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

import yaml

ROOT = Path(__file__).resolve().parents[1]


class ReleaseGuardTests(unittest.TestCase):
    def test_contact_ids_use_sorted_genomic_pixels(self):
        # Cooler is only needed during full matrix processing.
        src = ROOT / 'scripts' / 'build_contact50_counts.py'
        spec = importlib.util.spec_from_file_location('microc_contact_counts', src)
        mod = importlib.util.module_from_spec(spec)
        self.assertIsNotNone(spec.loader)
        with patch.dict(sys.modules, {'cooler': types.ModuleType('cooler')}):
            spec.loader.exec_module(mod)
        first = {'W1': {(7, 12): 5, (2, 5): 3}, 'K1': {(2, 5): 7, (6, 9): 1}}
        permuted = {'K1': {(6, 9): 1, (2, 5): 7}, 'W1': {(2, 5): 3, (7, 12): 5}}
        wanted = [(2, 5), (6, 9), (7, 12)]
        self.assertEqual(mod.ordered_candidate_keys(first), wanted)
        self.assertEqual(mod.ordered_candidate_keys(permuted), wanted)

    def test_methods_are_general_and_do_not_claim_unrun_sensitivity(self):
        with tempfile.TemporaryDirectory() as td:
            tmp_path = Path(td)
            config = yaml.safe_load((ROOT / 'config' / 'config.template.yaml').read_text())

            # Five sample libraries: exercise a non-fixed sample count.
            samples = {}
            for i in range(1, 4):
                samples[f'CONTROL_{i}'] = {'condition': 'CTRL', 'clone': f'CTRL_{i}'}
            for i in range(1, 3):
                samples[f'TREATMENT_{i}'] = {'condition': 'DRUG', 'clone': f'TX_{i}'}
            for name, info in samples.items():
                for read in ('r1', 'r2'):
                    fastq = tmp_path / f'{name}_{read}.fastq.gz'
                    with gzip.open(fastq, 'wt') as fh:
                        fh.write('@synthetic\n' + 'A' * 71 + '\n+\n' + 'I' * 71 + '\n')
                    info[read] = str(fastq)

            config['samples'] = samples
            config['comparison'] = {'control_condition': 'CTRL', 'test_condition': 'DRUG'}
            config['loops']['sensitivity_resolutions'] = [5000, 10000, 25000]
            cfg = tmp_path / 'test.yaml'
            cfg.write_text(yaml.safe_dump(config))
            md = tmp_path / 'methods.md'
            paper = tmp_path / 'methods.txt'
            subprocess.run([
                sys.executable, str(ROOT / 'scripts' / 'build_microc_methods.py'),
                '--config', str(cfg), '--pipeline-dir', str(ROOT),
                '--output-md', str(md), '--output-txt', str(paper),
            ], check=True, capture_output=True, text=True)

            detailed = md.read_text()
            manuscript = paper.read_text()
            self.assertIn('5 libraries', detailed)
            self.assertIn('3 CTRL', detailed)
            self.assertIn('2 DRUG', detailed)
            self.assertIn('all 5 samples', detailed)
            self.assertIn('every test sample', detailed)
            self.assertNotIn('GRCh38 primary assembly', manuscript)
            self.assertIn('Configured optional loop resolutions (execution not verified)', detailed)
            self.assertIn('5 kb and 25 kb', detailed)
            self.assertNotIn('sensitivity analysis was performed', manuscript)
            self.assertNotIn('Resolution sensitivity analysis', manuscript)


if __name__ == '__main__':
    unittest.main()
