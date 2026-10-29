# SnakeMicroC validation gates

The workflow is research software. The status of individual tests must be
recorded for each release; do not infer a production guarantee from synthetic
test success or from one representative experiment.

| Gate | Evidence required |
| --- | --- |
| Source | Canonical `rules/Snakefile`, working root symlink, generic config, no experiment data. |
| Python | Unit tests for config, CLI, reports, shallow mode, matrix reuse and QC views. |
| R/statistics | Exercise edgeR thresholds, normalization, low counts and modeling assumptions on fixtures. |
| Snakemake | DAG parses and dry-runs for QC/full modes with one and multiple samples. |
| Environments | Fresh Conda solve, documented versions and file compatibility. |
| Real data | Complete representative full-depth and shallow workflows with audited outputs. |
| Scientific QA | Mapping/dedup, cis/trans, complexity, HiCRep, structure tracks, distance filters, loop diagnostics. |
| Reporting | HTML links, index, per-sample counts, QC labels, Methods and provenance verified. |
| Release | Reviewed git diff, no paths/data/secrets, tested revision and no force push. |

## Local development checks

```bash
python3 -m unittest discover -s tests -p 'test*.py' -v
python3 -m compileall -q scripts bin
bin/microc validate --config /path/to/workspace/config.yaml --check-paths
bin/microc qc --config /path/to/workspace/config.yaml --workdir /path/to/workspace --cores 8 --dry-run
bin/microc run --config /path/to/workspace/config.yaml --workdir /path/to/workspace --cores 16 --dry-run
```

Synthetic Python tests do not exercise external aligners, peak callers,
Conda solving, R model fitting, or end-to-end processing of actual FASTQs.
For different genomes, independent adaptation and validation of reference
handling are required.

## Scientific cautions

- Group-level structural differences are descriptive unless formally tested.
- Clones, libraries and sequencing runs are not automatically independent
  biological replicates.
- A differential contact is not, on its own, proof of a specific enhancer-
  promoter interaction or formation/loss of a loop.
- Reusing existing matrices requires verified bin/pixel equivalence and
  provenance; never bypass the checks to avoid computational cost.
