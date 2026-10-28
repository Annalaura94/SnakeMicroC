# SnakeMicroC - operating guide for new experiments and additional sequencing

This guide supplements the root [README](../README.md). It describes the
currently implemented command-line interface; it does not imply that upstream
FASTQ aggregation or cross-library integration is automated.

## 1. Start a completely new experiment

Keep one reusable code checkout and one workspace **per analysis unit**. Never
place raw data under the Git repository. A new sequencing run should first be
recorded with the biological sample ID, library-preparation ID, sequencing-run
ID, FASTQ paths and assembly/reference versions.

```bash
SRC=/home/reference_data/pipelines/SnakeMicroC
WORK=/absolute/path/to/new/experiment
SNAKE=/absolute/path/to/working/snakemake
CONDA=/absolute/path/to/snakemake_conda_cache
mkdir -p "$WORK"
cp "$SRC/config/config.template.yaml" "$WORK/config.yaml"
# Edit all sample/R1/R2, reference, scratch and comparison fields.
"$SRC/bin/microc" validate --config "$WORK/config.yaml" --check-paths
"$SRC/bin/microc" run --config "$WORK/config.yaml" \
  --workdir "$WORK" --snakemake "$SNAKE" --conda-prefix "$CONDA" \
  --cores 16 --dry-run
# After reviewing the sample-to-FASTQ mapping and the planned DAG:
"$SRC/bin/microc" run --config "$WORK/config.yaml" \
  --workdir "$WORK" --snakemake "$SNAKE" --conda-prefix "$CONDA" \
  --cores 16
```

`run` is the full workflow. It can skip or refuse differential analyses depending
on `analysis.differential`, `analysis.require_differential` and the actual
replication in the config. The current implementation and references are
human/GRCh38-oriented; do not claim organism-independent validation.

For a **new pilot** (QC only), set `analysis.mode: shallow_qc`, validate the
config and call `microc qc`; the shallow command builds the HTML report and
Methods too. `scripts/init_shallow_qc.py` can initialize a new workspace from
an existing validated YAML and a folder of FASTQ files. Always inspect the
generated samples list: **the entire FASTQ directory may contain unrelated
sequencing libraries**.

## 2. Terminology: biological sample, library and sequencing run

These entities must not be conflated in the experiment manifest:

- Biological sample: original biological preparation/replicate. Keep its real
  experimental relationship and condition; different clones or preparations
  do not automatically constitute biological replication.
- Library: independently prepared Micro-C library from the relevant material.
  A library preparation can differ technically while using the same original
  biological sample or DNA material.
- Sequencing run: an acquisition of reads from a library. Sequencing the same
  library again yields additional reads, **not** a new library or biological
  replicate.

For example, a second sequencing run of an existing library set is a
resequencing event rather than a new replicate. Repreparing libraries from
existing experimental DNA is a distinct library preparation but does not
establish new biological replication. Always record run IDs and library IDs in
a project-local manifest rather than inferring them from directory names.

## 3. What to do when the sequencing provider sends more reads

### Additional reads from the **same library**

1. Preserve every original FASTQ and sequence/run/flowcell identifier;
   register the new FASTQs and the exact library ID.
2. Confirm read pairing, naming, assembly, mapping and filtering consistency
   before combining sources.
3. Design a new versioned integration workspace with all the correct inputs.
   If combining FASTQs, do so only with an explicitly validated aggregation
   plan, keeping the original unmodified FASTQ files as sources.
4. **Deduplicate across the aggregate sequencing observations.** Do **not**
   concatenate the `*.nodups.pairs.gz` outputs from two separate runs and
   treat them as a globally deduplicated library: overlapping molecules may
   be counted twice.
5. Recompute dependent matrices/QC/analysis with traceable config and
   reproducibility evidence. Existing mapping results can be reused only when
   the intermediate file format and cross-run deduplication logic support it
   and this reuse has been validated.

The CLI currently expects one R1/R2 FASTQ path per configured sample. It does
**not** perform an automatic, provenance-aware aggregation of new sequencing
runs. Manually changing a finished workspace YAML in place is discouraged:
create a new analysis version and preserve the preceding outputs.

### Reads from a **newly prepared library**

Keep the library IDs separate and QC the preparations individually. Define
whether the material represents biological replicates, technical preparations,
or other dependencies before pooling or fitting differential models. Separate
library IDs are not interchangeable with independent biological replication.
Do not assign new preparations to the same sample key merely to force an
existing four-sample contrast to run.

## 4. Re-run, reporting and provenance

The wrapper is `bin/microc` with commands `validate`, `qc`, `run`, `report`.
`report` may build missing dependencies: it is *not* a read-only export.
Always dry-run after modifying config, code or workspace. On new experiments,
use the default full rerun triggers. The alternative `--rerun-triggers mtime`
and `--skip-path-checks` exist for **audited legacy migrations only**;
otherwise they may conceal changed inputs.

Every analysis archive should retain at least:

- the exact YAML used, source FASTQ manifest, pipeline Git commit and reference;
- `results/report/` including Methods and output index;
- `results/qc/`, `results/provenance/` and scientific results actually used;
- links/paths and integrity/backup information for FASTQ, pairs, and matrices;
- original versioned workspace and its Snakemake state until resume and
  reproducibility have been verified and backed up.

Never infer that an archived report proves a workspace can be deleted.

## 5. Release and cleanup

See [RELEASE_CHECKLIST](RELEASE_CHECKLIST.md) for tests and Git steps.
Cleanup must be preceded by an inventory and storage verification. Delete
nothing from a live workspace or legacy integration folder merely because
it appears redundant. The project-specific staging tool creates **copies of
small final outputs only**; it is not a cleanup/deletion tool.
