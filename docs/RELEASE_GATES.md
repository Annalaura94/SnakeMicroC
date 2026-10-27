# SnakeMicroC RC2 — gates to stable release

| Gate | Requirement | RC2 status |
| --- | --- | --- |
| Source tree | One canonical `rules/Snakefile`, root symlink entrypoint; no copied, independently maintained rule versions | Prepared in isolated source snapshot |
| Python unit tests | Config, CLI, summary report, Methods, ID determinism | **Passed locally** (8 synthetic tests) |
| R model tests | Confirm contact 50-kb filters, edgeR CLI thresholds, loop model, TMM and low-count behavior using small fixtures | Not yet tested |
| Snakemake DAG | Parse from root symlink; dry-run `qc`, `report`, `all` with single/multi-sample configs | **Required on server** |
| Conda | Fresh resolution, record exact versions, compatibility across five env files | **Required on server** |
| End-to-end | Real MIAT run or audited reuse of existing outputs; check expected counts and file integrity | **Required on server** |
| Scientific QA | QC depth/complexity, HiCRep, compartment phase, insulation, contact-distance filtering, loop sensitivity | Not yet closed |
| Report | Correct HTML links and fully accurate Methods and per-sample metrics with real outputs | Synthetic tests passed; real MIAT still required |
| GitHub | Branch, reviewed diff, tests, no data/secrets, CI, protected merge/tag | Not published |
| Cleanup | Verified backup, hash ledger, retention policy, quarantine before deletion | Not authorized / not done |

## Checks to run on server

1. `python -m unittest discover -s tests -p 'test*.py' -v`
2. `bin/microc validate --config /path/to/MIAT/config_mustache_v1.yaml`
3. `bin/microc qc --config /path/to/config.yaml --workdir /path/to/experiment --cores 8 --dry-run --rerun-triggers mtime`
4. `bin/microc run --config /path/to/config.yaml --workdir /path/to/experiment --cores 24 --dry-run --rerun-triggers mtime`
5. Inspect the Snakemake DAG, missing files and rule triggers; record dry-run log in project `logs/`.
6. Only after review run targeted jobs with frozen code and documented parameters. Generate real report and compare it to frozen outputs.
7. When the full run is complete, repeat the dry-run to check no unexpected jobs remain.

## What is intentionally excluded from generic Micro-C core

MIAT-specific RNA-seq, ChIRP, triplex, cCRE and gene-locus interpretation belong in a separately versioned downstream integration module. They must not become hidden requirements for running a new Micro-C dataset.

## Reproducibility caveats

This snapshot is not an environment lock. `snakemake` was unavailable in the local container, and downloading it failed, so actual Snakemake parsing could not be certified here. Do not claim RC2 production-ready until the server checks above pass.


## RC2 matrix-reuse gate (27 Oct 2026)

- [x] Four MIAT full 50-kb raw pixel/bin comparisons PASS (reported by user).
- [ ] Verify `run --dry-run` on shared source using the audited MIAT manifest; no legacy `.mcool` rebuild may appear.
- [ ] Inspect all configured `.mcool` zoom levels and balancing weights on the server.
- [ ] Confirm single `insulation_cool` is not needed when insulation/loop resolutions agree.
- [ ] Finish requested 10-kb insulation, 50-kb differential contacts, QC and integrated report.
- [ ] Check that a second dry-run is no-op and that Methods accurately describe adopted matrices.
- [ ] Review Git diff; no backup directories, datasets or credentials in the commit.
