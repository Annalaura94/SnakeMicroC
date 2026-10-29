# SnakeMicroC

**SnakeMicroC** is a Snakemake workflow for paired-end Micro-C sequencing data.
It produces deduplicated chromatin contacts, QC summaries, contact matrices,
structural tracks, and, for appropriate sample designs, differential contacts
and candidate loops. An HTML report and a Methods summary are generated with
run-level provenance.

**Status:** Research workflow, under active validation. Tested on a human
reference workflow; new environments, designs, and non-human references need
independent validation. This repository contains pipeline source, not data or
experiment-specific analysis results.

## Quick start

Requirements: Linux, Python 3, a working Snakemake installation, and Conda for
rule-specific environments. Snakemake 7-style invocation was used during
validation; compatibility with other major versions has not been fully tested.
See [`env/`](env/) for package specifications. A fresh Conda solve is not a lockfile.

```bash
git clone https://github.com/Annalaura94/SnakeMicroC.git
cd SnakeMicroC
PIPELINE="$(pwd)"
WORK="$HOME/microc_workspaces/example"
mkdir -p "$WORK"
cp "$PIPELINE/config/config.template.yaml" "$WORK/config.yaml"
# Edit reference FASTA/FAI, tmpdir, samples (R1/R2), and group labels.

"$PIPELINE/bin/microc" validate --config "$WORK/config.yaml" --check-paths
"$PIPELINE/bin/microc" run --config "$WORK/config.yaml" \
  --workdir "$WORK" --cores 16 --dry-run
# After reviewing the sample mapping and dry-run:
"$PIPELINE/bin/microc" run --config "$WORK/config.yaml" \
  --workdir "$WORK" --cores 16
```

The wrapper assumes `snakemake` is already available on `PATH`. If it is not,
activate its environment or pass `--snakemake /path/to/snakemake`. Use
`--conda-prefix /path/to/cache` to share rule environments. For Snakemake
versions without `--conda-frontend`, add `--conda-frontend none`; actual
compatibility still needs a dry-run on that installation.

The example config describes a **two-control / two-test human** comparison.
It contains placeholder paths and **must** be edited before execution. Keep
config files, FASTQs, results, scratch files and logs in a workspace *outside*
the cloned source directory.

## Commands

| Command | Purpose |
| --- | --- |
| `bin/microc validate` | Validate config and, with `--check-paths`, reference and FASTQ paths. |
| `bin/microc qc` | Run QC targets; in `shallow_qc` mode also produce the QC report. |
| `bin/microc run` | Run the configured processing and analysis workflow. |
| `bin/microc report` | Build the report **and any missing dependencies**; not read-only. |

All commands accept `--config`; processing commands also accept `--workdir`,
`--cores`, `--dry-run`, `--snakemake` and `--conda-prefix`.

## Analysis modes and major outputs

**Full workflow:** BWA-MEM + pairtools processing, pairs QC, `.cool` / `.mcool`
matrices, P(s), reproducibility metrics, A/B compartments, insulation,
condition-blind loop candidates, and optional edgeR comparisons if the design
supports them. Set `analysis.differential: false` for a descriptive run;
`analysis.require_differential: true` rejects insufficient replication when
differential analysis is requested.

**Shallow pilot:** set `analysis.mode: shallow_qc` and
`analysis.differential: false`. Run `bin/microc qc` to generate mapping,
deduplication, complexity-screening and raw 1-Mb maps, report and Methods.
This is not a substitute for full-depth structural analysis. A helper can
prepare a fresh shallow config from the template and a FASTQ folder:

```bash
python "$PIPELINE/scripts/init_shallow_qc.py" \
  --base-config "$WORK/config.yaml" \
  --fastq-dir /path/to/paired_fastq_folder \
  --workspace "$HOME/microc_workspaces/pilot" \
  --condition PILOT
```

Check generated R1/R2 pair assignments before running. Multiple lanes are
rejected rather than silently merged. `PILOT` is only a generic group label.

Outputs are written under `<workdir>/results/`, including `pairs/`, `qc/`,
`mcool/`, `reproducibility/`, `structural/`, `contacts/`, `loops_mustache/`,
`provenance/` and `report/` as applicable to the chosen mode. The report
contains `MicroC_report.html`, `MicroC_summary.tsv` and Methods files.

## Important limitations

- Currently validated primarily with human chromosomes (`chr1`-`chr22`, `chrX`).
  Other references require code inspection and independent tests.
- Differential group statistics require an appropriate experimental design:
  technical library preparations are not automatically independent biological
  replicates. Descriptive structure comparisons are not statistical tests.
- The CLI accepts one FASTQ R1/R2 pair per configured sample; it does not
  automatically aggregate repeated sequencing runs or perform global dedup
  across independently processed read sets.
- Adapter trimming is not explicitly performed; BWA-MEM may soft-clip read
  ends, which is not an adapter-specific trimming step.
- Optional reuse of previously calculated `.mcool` files requires a verified
  equivalence report; never enable it for new data without an audit.

For details: [operating guide](docs/OPERATING_GUIDE.md),
[release checklist](docs/RELEASE_CHECKLIST.md),
[validation gates](docs/RELEASE_GATES.md) and [changelog](CHANGELOG.md).

## Developer checks

```bash
python3 -m unittest discover -s tests -p 'test*.py' -v
python3 -m compileall -q scripts bin
```

Synthetic tests do not replace end-to-end execution on real FASTQs, Snakemake
DAG parsing, Conda solving or statistical-model verification. The root
`Snakefile` is a symlink to `rules/Snakefile`.
