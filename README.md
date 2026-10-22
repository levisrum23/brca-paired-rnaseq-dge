# Breast cancer differential gene expression

**Question:** which genes and pathways differ between tumor and normal breast tissue in a public RNA-seq dataset?

**Data:** GEO [GSE133998](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133998): 6 patients, each with one tumor sample and one adjacent-normal sample (12 samples, raw read counts). Details, checksum and caveats are in [data/DATASET.md](data/DATASET.md).

**Result summary:** 2,823 genes differ (1,041 higher and 1,782 lower in tumor). Genes higher in tumor point to immune, cell-cycle and tissue-remodeling pathways; genes lower in tumor point to fat-cell and muscle-related pathways. Read the 2-page write-up in [report.md](report.md). With only 6 pairs, treat the results as exploratory.

## Rerun everything

Python 3.11 is needed. Create an environment from the pinned versions in `environment.yml` (conda), or with pip:

```bash
python3.11 -m venv .venv
.venv/bin/pip install pytest==9.1.1 pandas==3.0.6 numpy==2.4.6 scikit-learn==1.9.1 matplotlib==3.11.2 seaborn==0.13.2 pydeseq2==0.5.4 gseapy==1.3.1
```

Run the whole pipeline (downloads the 6 MB file the first time, then writes everything to `results/`):

```bash
.venv/bin/python run_all.py
```

Run the tests:

```bash
.venv/bin/pytest
```

The enrichment step needs internet (it fetches the gene-set libraries from Enrichr). Everything else works offline after the first download.

## Method

| Step | Script | What it does |
|---|---|---|
| 1 | `src/s01_download.py` | Download, SHA-256 check, load counts and metadata |
| 2 | `src/s02_qc.py` | Filter low-count genes, library sizes, PCA, correlation heatmap |
| 3 | `src/s03_dge.py` | PyDESeq2, design `~ patient + condition` (paired), tumor vs normal |
| 4 | `src/s04_plots.py` | Volcano, MA plot, top-gene heatmap |
| 5 | `src/s05_enrichment.py` | Over-representation test (Hallmark and KEGG gene sets from Enrichr), up and down genes separately, background = the genes tested |
| 6 | `src/s06_interpret.py` | Compare with 11 published marker directions (`data/known_markers.csv`, each with its source) |

**Settings** (in `src/config.py`): significant means adjusted p < 0.05 and |log2 fold change| > 1. Random seed 42. Low-count filter: a gene needs at least 10 reads in total and at least 1 read in 3 or more samples.

**Checks built in:** `run_all.py` stops if tumor and normal do not separate in the PCA (silhouette must be above 0). It warns if no gene is significant. 40 tests cover the functions, and one test reruns the whole pipeline on a small synthetic dataset.

## Limits

- Only 6 patients. Power is limited.
- "Normal" is tissue next to the tumor, not tissue from healthy donors.
- Bulk tissue: differences can come from a different mix of cell types (for example immune cells or fat cells), not only from changes inside tumor cells.
- The significance cutoff and the loose low-count filter let some noisy low-expression genes through (about a quarter of the significant genes). No fold-change shrinkage was applied.
- MMP11 could not be tested (2 reads in total). See `data/DATASET.md`.

## Folder

```
data/        DATASET.md, metadata.csv, known_markers.csv (raw download is not committed)
src/         one script per step
tests/       pytest tests
results/     tables and figures written by run_all.py
report.md    2-page write-up
docs/        design spec and implementation plan
```

## Licence

MIT, see [LICENSE](LICENSE). Data is from GEO GSE133998; please cite the original authors.
