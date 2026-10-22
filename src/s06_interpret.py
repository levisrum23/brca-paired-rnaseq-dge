import pandas as pd

from src.config import DATA_DIR, RESULTS_DIR


def compare_markers(res: pd.DataFrame, markers: pd.DataFrame) -> pd.DataFrame:
    """Compare observed tumor-vs-normal direction with the expected direction.

    agrees is True/False by the sign of log2FoldChange; None if the gene is absent
    or has no result. Significance is not part of the call: padj is shown beside it.
    """
    rows = []
    for m in markers.itertuples(index=False):
        lfc = padj = None
        agrees = None
        if m.gene in res.index and pd.notna(res.loc[m.gene, "log2FoldChange"]):
            lfc = float(res.loc[m.gene, "log2FoldChange"])
            padj = res.loc[m.gene, "padj"]
            padj = None if pd.isna(padj) else float(padj)
            agrees = bool((lfc > 0) == (m.expected_direction == "up"))
        rows.append({"gene": m.gene, "expected_direction": m.expected_direction, "observed_lfc": lfc,
                     "padj": padj, "agrees": agrees, "source": m.source})
    return pd.DataFrame(rows, columns=["gene", "expected_direction", "observed_lfc", "padj", "agrees", "source"]
                        ).astype({"agrees": object})


def main() -> None:
    res = pd.read_csv(RESULTS_DIR / "dge_all.csv", index_col=0)
    markers = pd.read_csv(DATA_DIR / "known_markers.csv")
    out = compare_markers(res, markers)
    out.to_csv(RESULTS_DIR / "marker_comparison.csv", index=False)
    n_true = sum(a is True for a in out["agrees"])
    n_false = sum(a is False for a in out["agrees"])
    n_none = sum(a is None for a in out["agrees"])
    print(f"markers: {n_true} agree, {n_false} disagree, {n_none} not testable")
    print(out[["gene", "expected_direction", "observed_lfc", "padj", "agrees"]].to_string(index=False))


if __name__ == "__main__":
    main()
