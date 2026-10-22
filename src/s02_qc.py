from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

from src.config import RESULTS_DIR, SEED

COLORS = {"normal": "#0072B2", "tumor": "#D55E00"}


def filter_low_counts(counts: pd.DataFrame, min_total: int = 10, min_samples: int = 3) -> pd.DataFrame:
    """Keep genes with total count >= min_total and count >= 1 in at least min_samples samples."""
    keep = (counts.sum(axis=0) >= min_total) & ((counts >= 1).sum(axis=0) >= min_samples)
    return counts.loc[:, keep]


def _pcs(counts: pd.DataFrame) -> np.ndarray:
    return PCA(n_components=2, random_state=SEED).fit_transform(np.log1p(counts.to_numpy(dtype=float)))


def pca_separation(counts: pd.DataFrame, meta: pd.DataFrame) -> float:
    """Silhouette of tumor/normal labels on the first 2 PCs of log1p counts (-1 to 1)."""
    return float(silhouette_score(_pcs(counts), meta.loc[counts.index, "condition"]))


def run_qc(counts: pd.DataFrame, meta: pd.DataFrame, results_dir: Path = RESULTS_DIR) -> pd.DataFrame:
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    filtered = filter_low_counts(counts)

    sizes = counts.sum(axis=1).rename("library_size").to_frame()
    sizes["possible_outlier"] = sizes["library_size"] < sizes["library_size"].median() / 3
    sizes.to_csv(results_dir / "qc_library_sizes.csv")

    cond = meta.loc[filtered.index, "condition"]
    pcs = _pcs(filtered)
    fig, ax = plt.subplots(figsize=(6, 5))
    for c, color in COLORS.items():
        m = (cond == c).to_numpy()
        ax.scatter(pcs[m, 0], pcs[m, 1], c=color, label=c, s=50)
    for (x, y), name in zip(pcs, filtered.index):
        ax.annotate(name, (x, y), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PCA of log1p counts")
    ax.legend()
    fig.tight_layout()
    fig.savefig(results_dir / "qc_pca.png", dpi=150)
    plt.close(fig)

    corr = np.log1p(filtered).T.corr(method="spearman")
    g = sns.clustermap(corr, cmap="viridis", figsize=(7, 7),
                       row_colors=cond.map(COLORS), col_colors=cond.map(COLORS))
    g.fig.suptitle("Sample correlation (Spearman, log1p)")
    g.savefig(results_dir / "qc_correlation.png", dpi=150)
    plt.close(g.fig)
    return filtered


def main() -> None:
    from src.s01_download import COUNTS_PATH, META_PATH, load_dataset

    counts, meta = load_dataset(COUNTS_PATH, META_PATH)
    filtered = run_qc(counts, meta)
    print(f"genes kept: {filtered.shape[1]} of {counts.shape[1]}; PCA silhouette: {pca_separation(filtered, meta):.3f}")


if __name__ == "__main__":
    main()
