# Changelog

## Repository presentation cleanup - 2026-10-29

- Separated project-specific migration notes from reusable pipeline source.
- Replaced local server paths in public documentation with portable examples.
- Changed the shallow-QC helper's default condition to generic `PILOT`; the
  condition may be set explicitly with `--condition`.
- Retained processing, matrix, statistics and reporting rules unchanged.

## RC2 shallow pilot, report and diagnostics - 2026-10-27/28

- Added QC-only shallow mode and a reusable FASTQ-pair configuration helper.
- Added full-depth descriptive structural diagnostics and additional report
  outputs with appropriate interpretation caveats.
- Improved run/version documentation and operating instructions.

## RC2 matrix handling and scaling updates - 2026-10-27

- Added opt-in reuse of independently audited legacy `.mcool` files.
- Reused existing verified contact matrices where possible and improved
  provenance of matrix adoption.
- Added independent BED4 views for `pairtools scaling` and regression tests.
- Extended CLI rerun-trigger controls for deliberate legacy migrations.

## RC2 workflow consolidation - 2026-10-26

- Consolidated processing, QC, contact matrices and structural analysis into
  one Snakemake entry point.
- Added deterministic contact IDs and configuration-driven differential
  thresholds; historical defaults may differ from newly configured thresholds.
- Added a config validator, runner, generic reports and Methods generation.
- Added synthetic Python regression tests.
