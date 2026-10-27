# MIAT Micro-C: controlled migration to the unified workflow

**26 Oct 2026:** The audit snapshot includes completed processing, QC, compartments, HiCRep, Mustache loops at 10 kb and a loop edgeR report, but does **not** list `results/contacts/` differential-contact outputs in the canonical `workspaces/v1` folder. Existing results are not evidence that the new integrated run is complete.

The MIAT experiment has four samples: `WT_Rep1`, `WT_Rep2`, `MIAT_D4`, `MIAT_A5`. Genome: hg38/GRCh38. Existing 10-kb pooled loop analysis: **9,535 candidate loops**, **1** primary edgeR FDR-significant; sensitivity to TMM and window choice requires cautious interpretation. Do not describe the two candidate loci as confirmed MIAT-regulated enhancer–promoter interactions.

## Server paths from the audit

- Shared source: `/home/reference_data/pipelines/SnakeMicroC`
- MIAT workspace: `/mnt/nvme_raid0/epigen/Maldotti_neurolncRNA/MicroC/workspaces/v1`
- MIAT config: `<workspace>/config_mustache_v1.yaml`
- Scratch/Conda environments: `/mnt/nvme_raid0/tmp/atamburr_microc_*`

## Workflow before changing anything

1. Check server `.screen` sessions, user processes and `git status`.
2. Create an independent release staging directory in `$HOME`; never apply files from this archive directly to the shared source tree or `workspaces/v1`.
3. Compare source files against the 26 October snapshot; abort if they have changed unexpectedly.
4. Execute `python -m unittest discover -s tests -p 'test*.py' -v` from staging.
5. Validate the real MIAT YAML; verify that the new configuration does not unexpectedly change `contact_maps`, `boundaries`, contact filters, loop quantification or differential parameters.
6. Run **dry-run only** using the original workspace. Compare the planned jobs with the prior Snakefile; inspect missing/modified outputs. Do not run a full rebuild without review.

Example (after RC copied into a staging directory, with existing Snakemake activated):

```bash
RC="$HOME/SnakeMicroC_rc2_20261026"
WORK="/mnt/nvme_raid0/epigen/Maldotti_neurolncRNA/MicroC/workspaces/v1"
"$RC/bin/microc" validate --config "$WORK/config_mustache_v1.yaml"
"$RC/bin/microc" run \
  --config "$WORK/config_mustache_v1.yaml" \
  --workdir "$WORK" \
  --cores 24 --dry-run --rerun-triggers mtime --skip-path-checks \
  --conda-prefix /mnt/nvme_raid0/tmp/atamburr_microc_snakemake_conda
```

The new `all` target may schedule structural 10-kb outputs and 50-kb differential contacts that do not exist under the current naming scheme. Compare their provenance and parameters before executing.

## Scientific acceptance checks (nonnegotiable)

- Verify every filter in the **new** 50-kb differential contact table is applied before model fitting; confirm distances in [150 kb, 5 Mb] and no diagonal pixels.
- Compare the primary 10-kb Mustache results against frozen reference outputs (9,535 candidates, one FDR-significant loop), explaining any changed counts.
- Verify sample libraries, normalized contact depth, HiCRep, E1, insulation and boundary figure regeneration.
- Confirm Methods/summary match actually completed jobs; ensure unexecuted 5-kb sensitivity is never claimed.
- Run fresh full DAG and subsequent no-op (or explained minimal-work) DAG after acceptable first run.

## Safe cleanup

Do not delete raw FASTQ, pairs, `.mcool`, 10-kb `.cool`, result tables, version snapshots, environments or `.snakemake` state before the new release passes the gates. Preserve previous analyses outside the canonical report as frozen legacy evidence; no file is categorized as safe to delete solely because its directory name looks temporary. See [release gates](RELEASE_GATES.md).


## 27 Oct 2026: full 50-kb raw-pixel equivalence established

All four samples (`WT_Rep1`, `WT_Rep2`, `MIAT_D4`, `MIAT_A5`) passed the
whole-genome bin-coordinate and raw-pixel SHA256 comparisons between their
existing `.mcool::/resolutions/50000` and the October 10-kb `.cool` data
coarsened to 50 kb. The recorded result is
`results/provenance/mcool_50kb_pixel_equivalence.tsv` (all `PASS`, exit 0).
This supports retaining the old `.mcool` at 50 kb. It is **not** a proof of
identity at every additional resolution or of balancing weights.

The consolidation patch allows **opt-in** adoption of these externally verified
`.mcool` files for this migrated experiment. It requires the complete PASS TSV,
all four sample IDs, the original 10-kb coolers and older input modification
times than the TSV evidence. The baseline template defaults to fresh `.mcool`
generation for new datasets.

To keep a single MIAT experiment configuration after the verified patch is
installed, preserve a backup of `config_mustache_v1.yaml` and add these two
values under its existing `contact_maps` dictionary (do not change other
analytical settings):

```yaml
  reuse_verified_existing_mcools: true
  existing_mcool_equivalence_tsv: /mnt/nvme_raid0/epigen/Maldotti_neurolncRNA/MicroC/workspaces/v1/results/provenance/mcool_50kb_pixel_equivalence.tsv
```

Run `bin/microc ... --dry-run` from the shared source path and verify that
**neither** `make_mcool` nor `make_base_cool` nor `make_qc_cool` is scheduled.
The updated insulation rule reuses the existing 10-kb per-sample loop cooler
when the configured insulation resolution is also 10 kb.

The migration must still pass all biological QC, 50-kb differential-contacts
review and full Methods/report checks before GitHub release or any cleanup.
