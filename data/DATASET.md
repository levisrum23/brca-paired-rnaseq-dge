# Dataset

- **Series:** GEO GSE133998, "Transcriptomic sequencing of breast cancer" (Homo sapiens, platform GPL20795, Illumina). Submitted by Xin Dong, Lei Wei et al. (contributors as listed on GEO).
- **Page:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133998
- **File used:** `GSE133998_cancer_and_adjacent_normal.txt.gz` (6,287,837 bytes, about 6 MB). Download URL is in `src/s01_download.py`.
- **SHA-256:** `3983a2180e3631bfc7141c7fd0241e9ae12679444b4a7012df583ec4ec08f590`. GEO publishes no checksum, so this is the value recorded at first download on 2026-10-22. It detects later corruption, not a wrong source.
- **Licence:** none stated on the series page. Cite the GEO accession when using it.
- **Not used:** the blood file `GSE133998_breast_cancer_blood_counts.txt.gz` (blood is not part of the tumor vs normal question).

## Design
- 6 patients (H36, H37, H38, H39, H40, H42). Each has one adjacent normal tissue sample (`-adj`) and one cancer tissue sample (`-cancer`): 6 normal + 6 tumor, fully paired.
- So the DE model should use `patient` as a covariate (paired design).
- **Limitation:** n = 6 pairs gives limited power. Results are exploratory. The report must say this.
- "Normal" here is tissue adjacent to the tumor, not healthy-donor tissue. It may carry some tumor-related changes.

## File format
- Tab-separated, one row per gene, columns: `Geneid, Chr, Start, End, Strand, Length`, then 12 sample columns of raw read counts (all whole numbers, checked in code).
- 61,155 gene rows. The loader drops 2 rows with no values (`MUC20-OT1`, `ANK2`) and merges 2 duplicate gene rows by summing, leaving 61,151 genes.
- Sample labels come from the GEO sample titles and are stored in `data/metadata.csv` (columns `sample, condition, patient`).

## QC observations (from `results/qc_*`)
- No library-size outliers (smallest is about 17.4 M reads, largest about 26.5 M).
- PCA silhouette for tumor vs normal is 0.326: the groups separate mostly along PC1.
- Patient H40's tumor and adjacent-normal samples sit close together in the PCA and cluster together in the correlation heatmap. This may mean low tumor purity or contamination in H40. It is a caveat for the report; the paired design handles patient effects.

## Sample removals
None. No sample was removed.

## Open caveat: MMP11
MMP11 has 2 reads in total across the 12 samples, so the low-count filter removed it and it cannot be tested. Published work (Fu 2015, BMC Cancer) reports it strongly up in tumors. The cause is unknown (possibly the gene annotation used to make this file). Report it as untestable, not as a disagreement.
