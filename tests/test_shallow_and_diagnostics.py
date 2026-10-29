"""Synthetic tests for new shallow/structural workflows, no real sequencing data."""
from __future__ import annotations
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from init_shallow_qc import discover_pairs, make_config
from structural_diagnostics import strict_switches, cluster_boundaries
from validate_config import ConfigurationError, load_and_validate

class ShallowTests(unittest.TestCase):
    def test_config_generation_pairs_paths_and_guards(self):
        with tempfile.TemporaryDirectory() as td:
            t=Path(td); fastq=t/'fastqs'; fastq.mkdir()
            for mate in ('1','2'):
                (fastq/f'Library_A_S1_L001_R{mate}_001.fastq.gz').write_bytes(b'example')
            template=yaml.safe_load((ROOT/'config/config.template.yaml').read_text())
            base=t/'base.yaml'; base.write_text(yaml.safe_dump(template))
            c=make_config(base,fastq,t/'newspace')
            self.assertEqual(c['analysis']['mode'],'shallow_qc')
            self.assertFalse(c['analysis']['differential'])
            self.assertFalse(c['contact_maps']['reuse_verified_existing_mcools'])
            self.assertEqual(list(c['samples']),['Library_A'])
            self.assertEqual(c['samples']['Library_A']['condition'], 'PILOT')
            configured = make_config(base, fastq, t/'other', condition='CONDITION_A')
            self.assertEqual(configured['samples']['Library_A']['condition'], 'CONDITION_A')
            with self.assertRaisesRegex(ValueError, 'Invalid condition label'):
                make_config(base, fastq, t/'invalid', condition='invalid name')
            conf=t/'config.yaml'; conf.write_text(yaml.safe_dump(c))
            self.assertEqual(load_and_validate(conf)['analysis_mode'],'shallow_qc')
            c['analysis']['differential']=True
            conf.write_text(yaml.safe_dump(c))
            with self.assertRaisesRegex(ConfigurationError,'shallow_qc'):
                load_and_validate(conf)
            (fastq/'Library_A_L002_R1_001.fastq.gz').write_bytes(b'other lane')
            with self.assertRaisesRegex(ValueError,'Multiple R1'):
                discover_pairs(fastq)

    def test_shallow_report_and_methods_do_not_claim_deep(self):
        with tempfile.TemporaryDirectory() as td:
            t=Path(td)
            cfg=yaml.safe_load((ROOT/'config/config.template.yaml').read_text())
            cfg['analysis']={'mode':'shallow_qc','differential':False,'require_differential':False}
            cfg['comparison']={};cfg['samples']={'PILOT': {'condition':'PILOT','r1':'absent_R1.fq.gz','r2':'absent_R2.fq.gz'}}
            conf=t/'c.yaml';conf.write_text(yaml.safe_dump(cfg))
            result=t/'results';(result/'report').mkdir(parents=True)
            (result/'qc').mkdir()
            (result/'qc/qc_status.txt').write_text('PILOT\tSCREEN_PASS')
            pd.DataFrame([{'sample':'PILOT','condition':'PILOT','no_dup_read_pairs':10000,
                           'dovetail_deep_qc_pass':'NOT_EVALUATED','shallow_screening_status':'SCREEN_PASS'}]).to_csv(result/'qc/qc_summary.tsv',sep='\t',index=False)
            for rel in ('qc/library_complexity.tsv','qc/scaling_curves.tsv','qc/contact_map_correlation_1Mb.tsv'):
                pd.DataFrame([{'sample':'PILOT','value':1}]).to_csv(result/rel,sep='\t',index=False)
            cool=result/'qc/matrices/PILOT.1000000.cool';cool.parent.mkdir(parents=True);cool.write_bytes(b'fake')
            for name in ('qc/scaling_curves.pdf','qc/library_complexity.pdf'):
                (result/name).write_bytes(b'PDF')
            (result/'provenance').mkdir()
            (result/'provenance/software_versions.tsv').write_text('software\tversion\n')
            (result/'provenance/resolved_config.yaml').write_text(conf.read_text())
            cmd=[sys.executable,str(ROOT/'scripts/build_microc_methods.py'),'--qc-only','--config',str(conf),'--pipeline-dir',str(ROOT),
                 '--output-md',str(result/'report/MicroC_methods_provenance.md'),
                 '--output-txt',str(result/'report/MicroC_methods_paper.txt')]
            r=subprocess.run(cmd,capture_output=True,text=True)
            self.assertEqual(r.returncode,0,r.stderr)
            methods=(result/'report/MicroC_methods_paper.txt').read_text()
            self.assertIn('Shallow',methods)
            self.assertIn('No deep multi-resolution',methods)
            self.assertNotIn('Loop discovery was performed',methods)
            r=subprocess.run([sys.executable,str(ROOT/'scripts/build_microc_report.py'),'--config',str(conf),
                  '--results',str(result),'--report-dir',str(result/'report'),'--differential','no','--shallow-qc'],capture_output=True,text=True)
            self.assertEqual(r.returncode,0,r.stderr)
            html=(result/'report/MicroC_report.html').read_text()
            self.assertIn('SHALLOW PILOT',html)
            self.assertNotIn('HiCRep SCC</h3>',html)
            self.assertIn('shallow_raw_qc_cool',(result/'report/MicroC_summary.tsv').read_text())

class StructuralTests(unittest.TestCase):
    def test_strict_e1_requires_both_samples_each_condition(self):
        df=pd.DataFrame([
            dict(chrom='chr1',start=0,end=100,CTRL1=.12,CTRL2=.21,KO1=-.2,KO2=-.11),
            dict(chrom='chr1',start=100,end=200,CTRL1=-.1,CTRL2=-.11,KO1=.3,KO2=.2),
            dict(chrom='chr1',start=200,end=300,CTRL1=.3,CTRL2=.3,KO1=-.2,KO2=.04),
        ])
        sw=strict_switches(df,['CTRL1','CTRL2'],['KO1','KO2'],.05)
        self.assertEqual(set(sw['switch']),{'A_to_B','B_to_A'})
        self.assertEqual(len(sw),2)
    def test_boundary_clusters_are_descriptive_and_strict(self):
        x=[pd.DataFrame([{'chrom':'chr1','center':100,'sample':'KO1'}]),
           pd.DataFrame([{'chrom':'chr1','center':105,'sample':'KO2'}]),
           pd.DataFrame([{'chrom':'chr1','center':400,'sample':'CTRL1'}]),
           pd.DataFrame([{'chrom':'chr1','center':408,'sample':'CTRL2'}])]
        out=cluster_boundaries(x,['CTRL1','CTRL2','KO1','KO2'],['CTRL1','CTRL2'],['KO1','KO2'],20)
        self.assertEqual(set(out.event),{'gain_in_test','loss_in_test'})
    def test_rules_exclude_structure_mustache_in_shallow(self):
        source=(ROOT/'rules/Snakefile').read_text()
        self.assertIn('if SHALLOW_QC:\n    ANALYSIS_TARGETS = []',source)
        self.assertIn('rule shallow_qc_map:',source)
        self.assertIn('if DIFFERENTIAL_ENABLED:\n    rule structural_diagnostics:',source)
        self.assertIn('    rule focal_candidate_figures:',source)

if __name__=='__main__':
    unittest.main()
