# SnakeMicroC operating guide

This guide describes the generic workflow, project structure and safe use of
additional sequencing data. Each experiment has its own config and workspace,
outside the cloned source repository.

## 1. Initialize an experiment

```bash
PIPELINE=/path/to/SnakeMicroC
WORK=/path/to/workspace
mkdir -p "$WORK"
cp "$PIPELINE/config/config.template.yaml" "$WORK/config.yaml"
# Edit genome FASTA/FAI, tmpdir, sample paths, conditions and analysis choices.
"$PIPELINE/bin/microc" validate --config "$WORK/config.yaml" --check-paths
"$PIPELINE/bin/microc" run --config "$WORK/config.yaml" \
  --workdir "$WORK" --cores 16 --dry-run
# Only after reviewing the planned jobs:
"$PIPELINE/bin/microc" run --config "$WORK/config.yaml" \
  --workdir "$WORK" --cores 16
```

You can pass `--snakemake /path/to/snakemake` and
`--conda-prefix /path/to/conda-cache` as appropriate for your installation.

For a shallow pilot, set `analysis.mode: shallow_qc` and
`analysis.differential: false`; use `microc qc`. The helper
`scripts/init_shallow_qc.py --base-config ... --fastq-dir ... --workspace ...`
can discover one paired R1/R2 FASTQ set per sample, but does **not** infer
biological conditions. It uses a generic `PILOT` label by default, adjustable
with `--condition`. Inspect the generated sample manifest carefully.

## 2. Samples, libraries and sequencing runs

Track separately in an **experiment-local manifest**:

- biological sample and independent-replicate relationships;
- library preparation identity, index set and insert/QC information;
- sequencing run/lane identifier, FASTQ source and checksums;
- reference genome, pipeline revision, config and software environment.

Resequencing an existing library produces additional reads, not a new library.
An additional library preparation from the same biological specimen is not
necessarily an independent biological replicate. Design differential analyses
using biologically justified replicates, not sample-name conventions.

## 3. More reads for an existing library

1. Retain original FASTQ files and record the new run/lane as a new source.
2. Check matching library IDs, indices, pairing, references and read lengths.
3. Build a new versioned integration workspace rather than overwriting a
   completed analysis.
4. Deduplicate **across the combined observations**. Do not concatenate
   separately deduplicated `nodups.pairs.gz` files and call the result globally
   deduplicated; the same molecule can occur in multiple sequencing runs.
5. Recompute downstream matrices, QC and analyses unless reuse is independently
   validated and documented. This pipeline does not automate cross-run merge.

For a new library preparation, keep library IDs separate and run individual QC.
Do not designate technical re-preparations as independent biological replicates
without an explicit experimental rationale.

## 4. Re-running and provenance

`bin/microc` accepts `validate`, `qc`, `run` and `report`. The last three use
Snakemake. `report` may execute dependencies: it is **not** a read-only export.
Always dry-run after editing code, config or source FASTQs. The default full
rerun triggers track input, parameter, code and timestamp changes. Alternative
triggers and `--skip-path-checks` are intended only for independently audited
legacy migrations; they must not conceal modified source data.

Archive the exact config, FASTQ manifest, source commit, reference, QC,
provenance, summary and final result files. Preserve the original analysis
workspace until restartability, checksums and backups are verified.

## 5. Boundaries of the workflow

The documented reference implementation is focused on human chromosomes.
Shallow screening does not perform deep A/B, insulation, loop or differential
analysis. Optional existing `.mcool` adoption requires explicit equivalence
evidence; current checks do not certify every balance weight or resolution.
See [release gates](RELEASE_GATES.md) before interpreting or publishing results.
