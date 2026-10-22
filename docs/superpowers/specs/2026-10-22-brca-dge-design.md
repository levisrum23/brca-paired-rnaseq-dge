# Breast Cancer Differential Gene Expression: Design

Date: 2026-10-22

## Purpose
A portfolio piece for a master's application in Medical Technology. It shows a clear biological question, a reproducible pipeline and a short honest write-up.

## Question
Which genes and pathways differ between tumor and normal breast tissue in a public RNA-seq dataset?

## Constraints
- Dry-lab only, $0 spend.
- Slow bandwidth: total download under about 50 MB. Use a processed count table, not raw reads.
- Language: Python only.
- Public data only. No access requests.

## Dataset (to be verified in step 1)
One public processed count table with tumor/normal labels and sample metadata (candidates: GEO, UCSC Xena, recount3). No accession is chosen yet. Step 1 must check size, licence, sample count and whether pairing or batch info exists. If nothing fits under the size limit, stop and report in BLOCKED.md.

## Layout
```
brca-dge/
  README.md          question, results, how to rerun
  environment.yml    pinned packages
  data/              raw (not committed) + metadata
  src/               one script per step
  results/           tables and figures
  report.md          2-page write-up
  tests/             pipeline and sanity checks
```

## Pipeline (one script per step)
1. Download and verify: fetch counts and metadata, record source URL and checksum.
2. QC: filter low-count genes, library sizes, PCA, sample-correlation heatmap. Remove outliers only with a written reason.
3. Differential expression: PyDESeq2, tumor vs normal. Add a pairing or batch covariate if the metadata has one.
4. Visuals: volcano, MA plot, heatmap of top genes.
5. Enrichment: gseapy on up and down genes, Hallmark and KEGG sets.
6. Interpretation: compare top genes with published breast cancer markers. Report agreements and disagreements.

## Rules
- Fixed random seed.
- Thresholds written in the README (default padj < 0.05, |log2FC| > 1).
- The report states only results the code produced.
- No invented links, accessions or numbers.

## Checks
- The pipeline reruns end to end from a clean folder.
- Sanity: PCA separates tumor from normal. Known markers (for example ESR1, ERBB2) behave plausibly.

## Out of scope (optional extras)
Survival analysis, machine-learning classifier, Streamlit dashboard.

## Success criteria
- One command reruns everything and regenerates all figures and tables.
- README and 2-page report are readable by a non-expert reviewer.
- Every claim in the report traces to a file in `results/`.
