# SnakeMicroC

An experimental, versioned **Micro-C analysis workflow** based on Snakemake. It processes paired-end FASTQ data into deduplicated pairs, multiresolution contact matrices, QC, structural tracks, and (when supported by the design) differential contacts and chromatin loops. It generates a navigable HTML report, TSV summary/index and Methods text from a single run target.

> **Status: audited RC2 development branch (28 October 2026).** The MIAT batch1 deep workspace completed 21/21 Snakemake steps on 27 October 2026, and the separate MIAT batch2 shallow-QC workspace completed 25/25 steps (four libraries). Python unit tests and MIAT/synthetic dry-runs were also verified during installation. This supports the tested MIAT workflows, **not** all possible genomes, sample designs or new experiments. Publication/merging of this branch remains a separate release step.

## Quick start

Create a dedicated project workspace. The repository contains source code only; keep FASTQ, BAM, pairs, coolers and results outside the Git repository.

```bash
mkdir -p /path/to/experiment
cp /path/to/SnakeMicroC/config/config.template.yaml /path/to/experiment/config.yaml
# Edit all reference, sample and scratch paths in config.yaml

/path/to/SnakeMicroC/bin/microc validate \
  --config /path/to/experiment/config.yaml --check-paths

/path/to/SnakeMicroC/bin/microc run \
  --config /path/to/experiment/config.yaml \
  --workdir /path/to/experiment \
  --cores 24 --dry-run

/path/to/SnakeMicroC/bin/microc run \
  --config /path/to/experiment/config.yaml \
  --workdir /path/to/experiment \
  --cores 24 --conda-prefix /path/to/snakemake-conda-cache
```

The `microc` command requires a separate Snakemake executable already on `PATH`. It does **not** install Snakemake or run in the background. On the laboratory system with an older Snakemake command use `--snakemake /absolute/path/to/snakemake` if necessary. The wrapper defaults to Snakemake 7-style `--use-conda --conda-frontend conda`, as used in the MIAT runs; use `--conda-frontend none` on newer versions that do not accept the frontend flag (actual compatibility requires testing).

| Command | Purpose |
| --- | --- |
| `microc validate` | Validate sample/configuration consistency; use `--check-paths` to verify reference and FASTQs |
| `microc qc` | Produce/readjust mapping, pairtools, complexity, contact-decay and map QC |
| `microc run` | Run all applicable processing, QC, structure, statistical tests, provenance and report |
| `microc report` | Build report and any **missing dependencies** (not a read-only export) |

For migration of already-produced MIAT results, temporarily use `--rerun-triggers mtime` and (only if FASTQ files have been archived) `--skip-path-checks` to avoid forcing a huge recalculation while comparing the old and new workflow. New experiments default to tracking file modification time, input, code and parameter changes. See [MIAT migration](docs/MIAT_MIGRATION.md) before running on MIAT.

## Analysis selection and limitations

- At least one sample: processing, QC and per-sample compartment/insulation tracks.
- At least two samples: HiCRep reproducibility.
- At least two control **and** two test samples: descriptive clone/group structural comparisons, edgeR 50-kb differential contacts and 10-kb differential loops (or configured resolutions).
- Other designs: differential results are marked **not run**, not falsely interpreted as absence of an effect. Set `analysis.require_differential: true` to stop when the design is insufficient.
- Genomic views are currently **human autosomes 1–22 and chrX**. Other organisms require adapting and validating `prepare_reference.py`; do not merely change `reference.assembly`.
- The QC contact maps/report figure names are fixed to the 1-Mb QC resolution. Keep `qc.qc_map_resolution: 1000000` in this RC.
- A/B E1 and insulation group-direction tables are **descriptive**, not statistical tests. A differential contact is not automatically a focal loop or evidence of enhancer–promoter causality.
- MIAT KO D4/A5 are distinct clones with one observation each; statistical significance depends on modeling assumptions and should not be described as validation across independent biological replicates.

## Pipeline stages

1. **Processing:** BWA-MEM, pairtools parse/sort/dedup, pairs stats.
2. **QC:** mapped reads, duplicates, cis/trans, library complexity, P(s), 1-Mb maps, correlation.
3. **Matrices and structural analysis:** cool/mcool, HiCRep, phased A/B compartment E1, cooltools insulation and boundary scores.
4. **Differential contacts:** autosomal cis contact candidates filtered by genomic distance **before** edgeR; parameterized count thresholds and FDR.
5. **Differential loops:** condition-blind pooled Mustache candidates, raw per-sample contact windows, edgeR quasi-likelihood, descriptive concordance.
6. **Report and provenance:** paper Methods, technical Methods, metric summaries, output index, version tables and resolved YAML.

## Main outputs

```text
experiment/results/
  pairs/                # per-sample deduplicated pairs/stats
  qc/                   # QC, complexity, scaling, map figures
  mcool/                # full-depth multiresolution maps
  reproducibility/      # HiCRep SCC
  structural/           # E1 and insulation tracks + group comparisons
  contacts/             # differential contacts (if eligible)
  loops_mustache/       # loop calls and differential results (if eligible)
  provenance/           # resolved config and software versions
  report/
    MicroC_report.html
    MicroC_summary.tsv
    MicroC_output_index.tsv
    MicroC_methods_paper.txt
    MicroC_methods_provenance.md
```

**Do not delete intermediates on the strength of a successful report alone.** A successful report establishes only that its expected input files exist and that corresponding Snakemake rules completed; it is not an independent assessment of scientific correctness, library quality or file integrity.

## Developer checks

```bash
python -m unittest discover -s tests -p 'test*.py' -v
python -m compileall -q scripts bin
```

The test suite uses fake Snakemake and synthetic TSVs; it does **not** test mapping, peak calling, R/edgeR execution, Conda solving or real Snakemake DAG compilation. See [release gates](docs/RELEASE_GATES.md).

## Reproducibility

- The Snakefile at repository root is a symlink to `rules/Snakefile`, not a duplicate entrypoint.
- `env/` contains reproducibility *specifications*, not lockfiles. Several dependencies remain range/unpinned; record exact versions for each run and lock tested environments before a final release.
- Every code change must be reviewed in Git with tests. Avoid project-specific scripts and exploratory outputs inside the generic pipeline.
- Example human config: `config/config.template.yaml`. Actual MIAT paths/config stay with the experiment, not in a public code repository.

For the scientific and migration caveats see [MIAT migration](docs/MIAT_MIGRATION.md) and [release gates](docs/RELEASE_GATES.md). For a new dataset or additional sequencing of an existing library, follow the [operating guide](docs/OPERATING_GUIDE.md) and [release checklist](docs/RELEASE_CHECKLIST.md).


### Audited reuse of legacy multi-resolution matrices

New datasets always build their own `.mcool` files. For an older dataset whose
existing `.mcool` matrices have been independently compared with newly generated
primary `.cool` matrices, set `contact_maps.reuse_verified_existing_mcools: true`
and provide `existing_mcool_equivalence_tsv` (absolute path). The entire sample
set must have `bins_match=YES`, `pixels_match=YES`, `status=PASS`. A report older
than either matrix is refused. In this mode the workflow intentionally has no
producer for the adopted `.mcool` files; missing files cause an error.

The adoption evidence currently certifies raw pixels at **50 kb only**, not all
matrix resolutions or balance weights. Inspect those separately during release
validation. For MIAT-specific instructions see `docs/MIAT_MIGRATION.md`.

`microc ... --rerun-triggers full` (the default) checks timestamps plus input,
code, parameters and Conda environment provenance. For a carefully audited
legacy migration only, `code-and-params` deliberately excludes timestamps.
This option must never be used to ignore unverified changed source data.


## Shallow-sequencing pilot QC (batch2 or any new sequencing batch)

The optional `analysis.mode: shallow_qc` runs mapping, deduplication, pairtools
statistics, P(s), occupancy-model library-complexity screening, and **raw,
unbalanced 1-Mb** QC maps. It does **not** execute deep multi-resolution
matrix balancing, HiCRep, E1, insulation, Mustache or differential testing.
It generates `results/report/MicroC_report.html` and both Methods outputs.
A screening pass is **not** a deep library QC pass. Sparse shallow contact-map
correlations and long-range complexity extrapolations are descriptive only.

Create a **separate** workspace from an existing validated config and a FASTQ
directory with unambiguously matched R1/R2 reads. Duplicate mates (e.g.,
multiple lanes) and incomplete pairs are rejected instead of silently merged.

```bash
SRC=/home/reference_data/pipelines/SnakeMicroC
BATCH2=/mnt/nvme_raid0/epigen/Maldotti_neurolncRNA/MicroC/workspaces/batch2_shallow
python "$SRC/scripts/init_shallow_qc.py" \
  --base-config /mnt/nvme_raid0/epigen/Maldotti_neurolncRNA/MicroC/workspaces/v1/config_mustache_v1.yaml \
  --fastq-dir /path/to/batch2_FASTQs \
  --workspace "$BATCH2"
"$SRC/bin/microc" validate --config "$BATCH2/config_shallow_qc.yaml" --check-paths
"$SRC/bin/microc" qc --config "$BATCH2/config_shallow_qc.yaml" \
  --workdir "$BATCH2" --cores 16 --snakemake /home/atamburr/.conda/envs/3plex/bin/snakemake \
  --conda-prefix /mnt/nvme_raid0/tmp/atamburr_microc_snakemake_conda --dry-run
```

Remove `--dry-run` only after reviewing the sample-to-FASTQ mapping and DAG.
The QC-only command automatically includes its report and Methods; it does
not process batch1 data. A change in FASTQ source layout, read depth, or
sequencing run always requires a new workspace and explicit provenance.
Do not treat newly prepared libraries as new biological replicates by default,
and do not simply concatenate separately deduplicated pairs when adding
sequencing reads from the same library. See `docs/OPERATING_GUIDE.md`.

## Added full-depth report diagnostics

For adequately replicated group comparisons, the integrated report now adds
strict descriptive E1 A/B sign switches (all control samples versus all test
samples, each above the configured absolute E1 threshold); insulation boundary
presence/absence clusters using the configured tolerance; and normalized-CPM
figures for the top core differential contact and loop candidates. These are
**descriptive** diagnostics, not extra statistical differential tests. They
must not be used to infer formation/loss of a loop or domain without
independent corroboration. The report links the coordinate tables and figures.
