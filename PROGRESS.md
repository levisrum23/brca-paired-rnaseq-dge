# Progress

Last updated: 2026-10-22

## Done
- Design spec and plan: `docs/superpowers/specs/`, `docs/superpowers/plans/`.
- Tasks 1 to 8 built and tested (39 tests pass). `python run_all.py` rebuilds all results from the 6 MB GEO file.
- Report (`report.md`) and README written. Every number in the report was checked against `results/`.
- Real-data result: 38,362 genes tested, 2,823 significant (1,041 up, 1,782 down in tumor).

## Open items
- Whole-branch review done (fresh reviewer, Opus). Two Important findings fixed (see git log). Deferred minors are listed below.
- MMP11 has 2 reads in this file; cause unknown (`data/DATASET.md`).
- Only 6 patients: results are exploratory. A second, larger dataset would strengthen the project.

## Deferred minor review notes
- NaN-padj genes are absent from volcano/MA plots (noted in the report).
- No fold-change shrinkage; low-count filter is loose (noted in the report).
- PCA uses log1p of raw counts (not size-normalized) and axes show no % variance.
- Download has no timeout; loader does not validate `condition` values or cast `patient` to string.

## Optional extras (not asked for)
- Validate the top genes in a larger dataset (for example a TCGA-BRCA processed table).
- Survival analysis, ML classifier, Streamlit dashboard (out of scope in the spec).
- Publish to GitHub for the portfolio (needs your OK: outward-facing).

## How to resume
Read `README.md`, then `data/DATASET.md`. Rerun with `.venv/bin/python run_all.py`. Environment: Python 3.11 venv (see README).
