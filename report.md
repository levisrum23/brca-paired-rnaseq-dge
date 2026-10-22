# Tumor vs normal gene expression in breast cancer: a paired RNA-seq analysis

*Dry-lab analysis of public data (GEO GSE133998). Code and rerun instructions: [README.md](README.md).*

## Question
Which genes and biological pathways differ between breast tumor tissue and the normal tissue next to it?

## Data
GEO series GSE133998: 6 breast cancer patients, each with one tumor sample and one adjacent-normal sample (12 samples). The file contains raw read counts for 61,155 gene rows. After removing 2 empty rows and merging 2 duplicate rows, 61,151 genes remained. Libraries had about 17 to 26 million reads each, with no outliers. Details: [data/DATASET.md](data/DATASET.md).

## Methods
1. **Filtering.** Genes needed at least 10 reads in total and at least 1 read in 3 or more samples. This left 38,362 genes.
2. **Quality control.** PCA and a sample-correlation heatmap (figures `results/qc_pca.png`, `results/qc_correlation.png`).
3. **Differential expression.** PyDESeq2 with the design `~ patient + condition`, which compares each tumor with its own patient's normal tissue. Positive log2 fold change means higher in tumor. A gene is called significant at adjusted p < 0.05 and |log2 fold change| > 1.
4. **Pathways.** Enrichr over-representation test on the up genes and the down genes separately, with the MSigDB Hallmark 2020 and KEGG 2026 gene sets.
5. **Markers.** Direction of change compared with 12 genes reported in three published studies (`data/known_markers.csv`).

Random seed 42. The whole analysis reruns with `python run_all.py`.

## Results

**Quality control.** In the PCA, tumor and normal samples separate mostly along PC1 (silhouette score 0.326, where 0 would mean no separation). One patient, H40, is an exception: its tumor and normal samples sit close together. I kept it, and the paired design accounts for patient differences.

**Differential expression.** 38,362 genes were tested. 2,823 were significant: **1,041 higher and 1,782 lower in tumor** ([volcano plot](results/volcano.png), [MA plot](results/ma_plot.png)). The most significant genes higher in tumor include COMP, LRRC15, MATN3, GAPT and EPYC (adjusted p < 0.0001 each). The most significant lower genes include VLDLR-AS1, TAT, EGFR-AS1, AC005089.1 and HIF3A (adjusted p < 0.0001 each). The top 30 genes are shown in the [heatmap](results/top_genes_heatmap.png).

![Volcano plot](results/volcano.png)

**Pathways.** Genes higher in tumor were enriched for 17 Hallmark sets and 25 KEGG sets (adjusted p < 0.05). The strongest Hallmark sets are Allograft Rejection, Interferon Gamma Response, Interferon Alpha Response, G2-M Checkpoint and Epithelial Mesenchymal Transition. The strongest KEGG sets are Cell Adhesion Molecule Interaction, Primary Immunodeficiency and Cytokine-Cytokine Receptor Interaction. **Genes lower in tumor had no significant Hallmark or KEGG set** (0 of each).

**Comparison with published markers.** Of 12 markers, 11 could be tested. All 11 changed in the published direction (0 disagreed). 7 of the 11 were also significant here: COL10A1, TUBB3, MMP1, COL11A1, S100P, LRRC15 and LEP. The other four (ADIPOQ, DLK1, OXTR, PIGR) moved in the expected direction but were not significant (adjusted p between 0.068 and 0.334), so I do not count them as confirmed. MMP11 could not be tested: it has 2 reads in total in this file and was removed by the filter. Published work reports it strongly higher in tumors, and I do not know why it is nearly absent here.

## Interpretation
The strongest signals are a rise in cell-cycle genes (G2-M Checkpoint), tissue-remodeling genes (Epithelial Mesenchymal Transition, collagens such as COL10A1 and COL11A1) and immune genes (interferon response, cytokine signaling). The immune result may partly mean more immune cells in tumor tissue, not only changes inside tumor cells, because this is bulk tissue. In the same way, a fall in fat-cell genes such as LEP may mean less fat tissue in the tumor sample.

## Limitations
- **Small sample.** 6 patients means low power. Results are exploratory and should be checked in a larger dataset.
- **Control tissue.** "Normal" is tissue next to the tumor. It can carry some tumor-related changes.
- **Patient H40.** Its tumor and normal samples look alike, which may mean low tumor purity or contamination.
- **Enrichment background.** Enrichr tests against all genes in each gene set, not against the genes expressed here. P-values may be slightly optimistic.
- **Not tested.** No survival analysis, no tumor subtypes, no validation in a second dataset.

## Where each number comes from
| Number | File |
|---|---|
| Samples, reads per library, outliers | `results/qc_library_sizes.csv` |
| Silhouette 0.326, genes kept (38,362) | printed by `python run_all.py` and `python -m src.s02_qc`; PCA in `results/qc_pca.png` |
| 38,362 tested, 2,823 significant, 1,041 up, 1,782 down | `results/dge_all.csv`, `results/dge_significant.csv` |
| Enriched sets (17 Hallmark, 25 KEGG, 0 for down genes) | `results/enrichment_*.csv` |
| Marker agreement (11 of 11 tested, 7 significant) | `results/marker_comparison.csv` |

## Sources of the marker directions
- Fu J, et al. 2015. Gene expression profiling leads to discovery of correlation of matrix metalloproteinase 11 and heparanase 2 in breast cancer progression. BMC Cancer 15:473. doi:10.1186/s12885-015-1410-y
- Makoukji J, et al. 2016. Gene expression profiling of breast cancer in Lebanese women. Scientific Reports 6:36639. doi:10.1038/srep36639
- Zhang M, et al. 2020. Bioinformatics analysis of prognostic significance of COL10A1 in breast cancer. Bioscience Reports 40(2):BSR20193286 (PMC7029149).

Data: Dong X, Wei L, et al., GEO GSE133998 (contributors as listed on the GEO page).
