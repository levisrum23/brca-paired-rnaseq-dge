import numpy as np
import pandas as pd
import pytest

from src.config import LOG2FC_MIN, PADJ_MAX
from src.s04_plots import ma_plot, top_gene_heatmap, volcano


@pytest.fixture
def res(synthetic):
    counts, _ = synthetic
    rng = np.random.default_rng(0)
    n = counts.shape[1]
    return pd.DataFrame({
        "baseMean": rng.uniform(10, 1000, n),
        "log2FoldChange": rng.normal(0, 2, n),
        "pvalue": rng.uniform(0, 1, n),
        "padj": rng.uniform(0, 1, n),
    }, index=counts.columns)


@pytest.fixture
def res_none_significant(res):
    return res.assign(padj=0.9)


def _nonempty(p):
    return p.exists() and p.stat().st_size > 0


def test_volcano_writes_png(res, tmp_path):
    volcano(res, PADJ_MAX, LOG2FC_MIN, tmp_path / "v.png")
    assert _nonempty(tmp_path / "v.png")


def test_ma_plot_writes_png(res, tmp_path):
    ma_plot(res, tmp_path / "ma.png")
    assert _nonempty(tmp_path / "ma.png")


def test_heatmap_writes_png(synthetic, res, tmp_path):
    counts, meta = synthetic
    top_gene_heatmap(counts, meta, res, 10, tmp_path / "h.png")
    assert _nonempty(tmp_path / "h.png")


def test_plots_handle_zero_significant(synthetic, res_none_significant, tmp_path):
    counts, meta = synthetic
    volcano(res_none_significant, PADJ_MAX, LOG2FC_MIN, tmp_path / "v.png")
    ma_plot(res_none_significant, tmp_path / "ma.png")
    top_gene_heatmap(counts, meta, res_none_significant, 10, tmp_path / "h.png")
    assert all(_nonempty(tmp_path / f) for f in ("v.png", "ma.png", "h.png"))


def test_plots_tolerate_nan_padj(synthetic, res, tmp_path):
    counts, meta = synthetic
    res = res.copy()
    res.iloc[:5, res.columns.get_loc("padj")] = np.nan
    res.iloc[:5, res.columns.get_loc("pvalue")] = np.nan
    volcano(res, PADJ_MAX, LOG2FC_MIN, tmp_path / "v.png")
    top_gene_heatmap(counts, meta, res, 10, tmp_path / "h.png")
    assert _nonempty(tmp_path / "h.png")
