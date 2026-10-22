from pathlib import Path

import matplotlib
import matplotlib.patches

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from src.config import LOG2FC_MIN, PADJ_MAX, RESULTS_DIR

UP, DOWN, NS = "#D55E00", "#0072B2", "#999999"
GRID = {"color": "#dddddd", "linewidth": 0.6}


def _status(res: pd.DataFrame, padj_max: float, lfc_min: float) -> pd.Series:
    sig = (res["padj"] < padj_max) & (res["log2FoldChange"].abs() > lfc_min)
    return pd.Series(np.where(~sig, "ns", np.where(res["log2FoldChange"] > 0, "up", "down")), index=res.index)


def _scatter_by_status(ax, x, y, status: pd.Series) -> None:
    for name, color, label in (("ns", NS, "not significant"), ("down", DOWN, "lower in tumor"), ("up", UP, "higher in tumor")):
        m = (status == name).to_numpy()
        ax.scatter(x[m], y[m], s=6, c=color, alpha=0.6 if name == "ns" else 0.8, linewidths=0,
                   label=f"{label} ({int(m.sum())})")
    ax.grid(**GRID)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, markerscale=2, fontsize=8)


def volcano(res: pd.DataFrame, padj_max: float, lfc_min: float, out: Path) -> None:
    res = res.dropna(subset=["padj", "log2FoldChange"])
    fig, ax = plt.subplots(figsize=(6, 5))
    _scatter_by_status(ax, res["log2FoldChange"].to_numpy(), -np.log10(res["padj"].clip(lower=1e-300)).to_numpy(),
                       _status(res, padj_max, lfc_min))
    ax.axhline(-np.log10(padj_max), ls="--", c="#555555", lw=0.8)
    for x in (-lfc_min, lfc_min):
        ax.axvline(x, ls="--", c="#555555", lw=0.8)
    ax.set_xlabel("log2 fold change (tumor vs normal)")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title("Volcano plot")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def ma_plot(res: pd.DataFrame, out: Path) -> None:
    res = res.dropna(subset=["padj", "log2FoldChange", "baseMean"])
    fig, ax = plt.subplots(figsize=(6, 5))
    _scatter_by_status(ax, np.log10(res["baseMean"] + 1).to_numpy(), res["log2FoldChange"].to_numpy(),
                       _status(res, PADJ_MAX, LOG2FC_MIN))
    ax.axhline(0, c="#555555", lw=0.8)
    ax.set_xlabel("log10 mean normalized count")
    ax.set_ylabel("log2 fold change (tumor vs normal)")
    ax.set_title("MA plot")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def top_gene_heatmap(counts: pd.DataFrame, meta: pd.DataFrame, res: pd.DataFrame, n: int, out: Path) -> None:
    """Top n genes by padj. Library-size normalized, log1p, z-scored per gene."""
    top = res.dropna(subset=["padj"]).nsmallest(n, "padj").index
    norm = counts.div(counts.sum(axis=1), axis=0) * counts.sum(axis=1).median()
    z = np.log1p(norm[top])
    z = (z - z.mean()) / z.std(ddof=0).replace(0, 1)
    order = meta.loc[z.index].sort_values(["condition", "patient"] if "patient" in meta else ["condition"]).index
    colors = meta.loc[order, "condition"].map({"normal": DOWN, "tumor": UP})
    g = sns.clustermap(z.loc[order].T, col_cluster=False, cmap="RdBu_r", center=0, col_colors=colors,
                       figsize=(7, 0.25 * len(top) + 2.5), yticklabels=True, cbar_kws={"label": "z-score"})
    g.fig.suptitle(f"Top {len(top)} genes by adjusted p-value", y=1.01)
    g.fig.legend(handles=[matplotlib.patches.Patch(color=DOWN, label="normal"),
                          matplotlib.patches.Patch(color=UP, label="tumor")],
                 loc="upper left", bbox_to_anchor=(0.02, 0.8), frameon=False, title="sample")
    g.savefig(out, dpi=150)
    plt.close(g.fig)


def main() -> None:
    from src.s01_download import COUNTS_PATH, META_PATH, load_dataset
    from src.s02_qc import filter_low_counts

    counts, meta = load_dataset(COUNTS_PATH, META_PATH)
    counts = filter_low_counts(counts)
    res = pd.read_csv(RESULTS_DIR / "dge_all.csv", index_col=0)
    volcano(res, PADJ_MAX, LOG2FC_MIN, RESULTS_DIR / "volcano.png")
    ma_plot(res, RESULTS_DIR / "ma_plot.png")
    top_gene_heatmap(counts, meta, res, 30, RESULTS_DIR / "top_genes_heatmap.png")
    print("wrote volcano.png, ma_plot.png, top_genes_heatmap.png")


if __name__ == "__main__":
    main()
