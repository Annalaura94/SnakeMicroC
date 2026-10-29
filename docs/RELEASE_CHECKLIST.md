# SnakeMicroC release checklist

This checklist concerns the **pipeline source**, not an individual project's
results. Use a feature branch and a reviewed commit; do not force-push.

1. Record `git status -sb`, current commit and remote. Start with a clean tree.
2. Confirm `git ls-files` contains only source, environment specifications,
   tests, templates and generic documentation. Do not commit local configs,
   experiment identifiers, personal paths, raw reads, matrices or reports.
3. Run `python3 -m unittest discover -s tests -p 'test*.py' -v` and
   `python3 -m compileall -q scripts bin` in a suitable environment.
4. Validate a synthetic or independent configuration. Verify Snakemake
   `qc --dry-run` and `run --dry-run` from a fresh workspace.
5. Verify a representative full-depth end-to-end run and a distinct shallow
   QC run. Record exact software versions, inputs, outputs, Methods and the
   code revision. A pass on one experiment does not certify every design.
6. Confirm that a no-change rerun does not silently schedule heavy matrix,
   mapping, differential or loop calls.
7. Review `git diff --check`, the complete staged diff and test results.
8. Push the reviewed feature branch and merge via an ordinary PR. Avoid
   history rewriting or forced updates to `main`.

Example:

```bash
cd /path/to/SnakeMicroC
git status -sb
git ls-files
git diff --check
python3 -m unittest discover -s tests -p 'test*.py' -v
python3 -m compileall -q scripts bin
# After review, commit only the intended source changes.
```

A working `.gitignore` does not untrack files previously committed; inspect
`git ls-files`, not just `git status`. Data and archive cleanup are separate
operations and must not be coupled to source releases.
