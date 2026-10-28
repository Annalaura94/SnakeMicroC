# SnakeMicroC - release checklist

Before updating the shared GitHub repository:

1. Confirm the installed source directory `/home/reference_data/pipelines/SnakeMicroC`
   is the single active code checkout, not a data workspace. Record commit,
   branch, `git status` and configured Git remote. Do not add outputs, raw data,
   caches or local lab configs.
2. Run `python -m unittest discover -s tests -p 'test*.py' -v` and
   `python -m compileall -q scripts bin` in the validated Python environment.
3. Confirm at least one representative full-depth workflow completed with
   its expected report, structural outputs and provenance; verify successful
   exit status. Separately confirm a representative shallow QC-only workflow
   completed with exactly the intended libraries.
4. On an unrelated synthetic/new dataset, validate the sample/FASTQ mapping,
   and complete a Snakemake dry-run. An earlier test run on MIAT does not by
   itself establish arbitrary-design or nonhuman support.
5. Re-run a **dry-run** of the completed reference workspace. It should not
   schedule heavy mapping, matrix creation, Mustache or differential analyses
   without a deliberate change to inputs or provenance.
6. Review `git diff --check`, `git status --short`, the changed documentation
   and exact staged paths. Push only the tested, reviewed Git commit, and
   use a feature branch/PR rather than forcing changes to `main`.

Example publish steps (only after steps 1-6 pass):

```bash
SRC=/home/reference_data/pipelines/SnakeMicroC
cd "$SRC"
git status -sb
git branch --show-current
git remote -v
git log -3 --oneline
git diff --check
# Review changes and commit intentional source/doc files only.
git push -u origin feature/microc-unified-rc2
```

If the remote branch diverged, stop and review; do not force-push.
Push only updates the feature branch. A merge to the repository's main branch
and a release/tag require their own review.

The repository `.gitignore` should continue to exclude FASTQ, pair files,
matrices and run outputs; test `git ls-files` before publishing, not just
`git status`, because already tracked files are unaffected by `.gitignore`.

**Archive and cleanup:** stage results and record checksums, verify source
locations and backup, compare manifests, and only then separately choose
nonessential duplicates for deletion. This release procedure performs no
destructive cleanup.
