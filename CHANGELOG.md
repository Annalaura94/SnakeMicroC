# Changelog

## RC2 matrix reuse consolidation (27 Oct 2026)

- For audited legacy datasets, opt-in verification of the complete 50-kb bin/pixel equivalence report lets Snakemake treat existing `.mcool` files as external inputs rather than silently rebuild them.
- For new datasets, generate `.mcool` atomically in a single rule directly from pairs, avoiding missing temporary `.cool` dependencies.
- Reuse QC-resolution `.mcool` URIs rather than independently reloading 1-Mb contact matrices.
- Reuse the primary loop-resolution `.cool` for insulation when resolutions match.
- Fix CLI rerun-trigger choices: full by default; code-and-params excludes mtime, explicitly restricted to audited migrations.
- Add adoption provenance to Methods/report and four regression tests.
- Server DAG and new processing path still need on-cluster validation before release.


## RC2 candidate — 2026-10-26 (NOT YET RELEASED)

- Consolidated the server audit source and original Conda YAMLs in one isolated staging tree.
- Added deterministic 50-kb contact IDs via genomic pixel sorting.
- Connected edgeR contact count thresholds to configuration (rather than hard-coded R values).
- Removed MIAT-specific text and fixed 3×3 assumptions from the generic loop report.
- Added a config validator, a `bin/microc` runner and conditional execution for designs lacking sufficient replication.
- Added an integrated HTML report, summaries and canonical output index; Methods are downstream of completed analysis targets.
- Corrected Methods wording for sample number, condition labels and unexecuted sensitivity analyses.
- Added synthetic Python regression tests and migration/release gates.

**Methods change:** new configs with `contacts.edgeR_min_total_count: 15` will actually use 15, whereas the historical R script used 20 unconditionally. The MIAT config from the audit omits this key; the new Snakefile retains the historical fallback of 20. Any scientific re-analysis after changing this parameter requires a new run and documented comparison.

**No data processed, deleted or published by this RC2 package.** Real Snakemake/Conda/R validation is outstanding.
