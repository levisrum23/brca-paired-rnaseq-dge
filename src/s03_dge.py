import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from src.config import LOG2FC_MIN, PADJ_MAX, RESULTS_DIR


def run_dge(counts: pd.DataFrame, meta: pd.DataFrame, covariates: list[str] | None = None) -> pd.DataFrame:
    """PyDESeq2 tumor vs normal. Positive log2FoldChange means higher in tumor."""
    missing = [s for s in counts.index if s not in meta.index]
    if missing:
        raise ValueError(f"samples missing from metadata: {missing}")
    meta = meta.loc[counts.index]
    design = "~ " + " + ".join([*(covariates or []), "condition"])
    dds = DeseqDataSet(counts=counts, metadata=meta, design=design, n_cpus=1, quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["condition", "tumor", "normal"], n_cpus=1, quiet=True)
    stats.summary()
    return stats.results_df


def significant(res: pd.DataFrame, padj_max: float, lfc_min: float) -> pd.DataFrame:
    return res[(res["padj"] < padj_max) & (res["log2FoldChange"].abs() > lfc_min)]


def is_fully_paired(meta: pd.DataFrame) -> bool:
    """True if every patient has at least one tumor and one normal sample."""
    conds = meta.groupby("patient")["condition"].agg(set)
    return bool(conds.map(lambda c: {"tumor", "normal"} <= c).all())


def main() -> None:
    from src.s01_download import COUNTS_PATH, META_PATH, load_dataset
    from src.s02_qc import filter_low_counts

    counts, meta = load_dataset(COUNTS_PATH, META_PATH)
    counts = filter_low_counts(counts)
    covariates = ["patient"] if "patient" in meta.columns and is_fully_paired(meta) else []
    res = run_dge(counts, meta, covariates)
    sig = significant(res, PADJ_MAX, LOG2FC_MIN)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    res.to_csv(RESULTS_DIR / "dge_all.csv")
    sig.to_csv(RESULTS_DIR / "dge_significant.csv")
    print(f"design covariates: {covariates}; tested {len(res)} genes; significant: {len(sig)} "
          f"({int((sig['log2FoldChange'] > 0).sum())} up, {int((sig['log2FoldChange'] < 0).sum())} down)")


if __name__ == "__main__":
    main()
