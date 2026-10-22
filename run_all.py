"""Run the whole pipeline: python run_all.py"""
from pathlib import Path

import pandas as pd

from src.config import DATA_DIR, LOG2FC_MIN, PADJ_MAX, RESULTS_DIR
from src.s01_download import COUNTS_PATH, META_PATH, load_dataset, main as download
from src.s02_qc import pca_separation, run_qc
from src.s03_dge import is_fully_paired, run_dge, significant
from src.s04_plots import ma_plot, top_gene_heatmap, volcano
from src.s05_enrichment import LIBRARIES, run_enrichment, split_up_down
from src.s06_interpret import compare_markers


def check_sanity(pca_sep: float, n_significant: int) -> None:
    """Hard failure if tumor and normal do not separate; only a warning if nothing is significant."""
    if pca_sep <= 0:
        raise RuntimeError(f"PCA sanity check failed: tumor/normal silhouette is {pca_sep:.3f} (must be > 0)")
    if n_significant == 0:
        print("warning: no genes passed the significance thresholds")


def run_pipeline(counts_path: Path, meta_path: Path, markers_path: Path, results_dir: Path) -> None:
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    counts, meta = load_dataset(counts_path, meta_path)
    filtered = run_qc(counts, meta, results_dir=results_dir)
    sep = pca_separation(filtered, meta)

    covariates = ["patient"] if "patient" in meta.columns and is_fully_paired(meta) else []
    res = run_dge(filtered, meta, covariates)
    sig = significant(res, PADJ_MAX, LOG2FC_MIN)
    res.to_csv(results_dir / "dge_all.csv")
    sig.to_csv(results_dir / "dge_significant.csv")
    check_sanity(sep, len(sig))

    volcano(res, PADJ_MAX, LOG2FC_MIN, results_dir / "volcano.png")
    ma_plot(res, results_dir / "ma_plot.png")
    top_gene_heatmap(filtered, meta, res, 30, results_dir / "top_genes_heatmap.png")

    up, down = split_up_down(sig)
    for direction, genes in (("up", up), ("down", down)):
        for key, library in LIBRARIES.items():
            run_enrichment(genes, library).to_csv(results_dir / f"enrichment_{direction}_{key}.csv", index=False)

    compare_markers(res, pd.read_csv(markers_path)).to_csv(results_dir / "marker_comparison.csv", index=False)
    print(f"done: PCA silhouette {sep:.3f}; {len(res)} genes tested; {len(sig)} significant "
          f"({int((sig['log2FoldChange'] > 0).sum())} up, {int((sig['log2FoldChange'] < 0).sum())} down); "
          f"covariates {covariates}")


if __name__ == "__main__":
    download()
    run_pipeline(COUNTS_PATH, META_PATH, DATA_DIR / "known_markers.csv", RESULTS_DIR)
