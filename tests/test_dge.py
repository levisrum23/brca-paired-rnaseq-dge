import numpy as np
import pandas as pd
import pytest

from src.config import LOG2FC_MIN, PADJ_MAX, SEED
from src.s03_dge import run_dge, significant


@pytest.fixture(scope="module")
def planted():
    """12 paired samples x 60 genes. u0-u9 are 8x up in tumor, d0-d9 are 8x down."""
    rng = np.random.default_rng(SEED)
    samples = [f"P{i}-{t}" for i in range(6) for t in ("adj", "cancer")]
    tumor = np.array([s.endswith("cancer") for s in samples])
    mu = np.full((12, 60), 500.0)
    mu[np.ix_(tumor, range(10))] *= 8
    mu[np.ix_(tumor, range(10, 20))] /= 8
    counts = pd.DataFrame(rng.poisson(mu), index=samples,
                          columns=[f"u{i}" for i in range(10)] + [f"d{i}" for i in range(10)] + [f"n{i}" for i in range(40)])
    meta = pd.DataFrame({"condition": np.where(tumor, "tumor", "normal"),
                         "patient": [s.split("-")[0] for s in samples]}, index=samples)
    return counts, meta


@pytest.fixture(scope="module")
def result(planted):
    counts, meta = planted
    return run_dge(counts, meta, covariates=["patient"])


def test_planted_up_genes_have_positive_lfc_and_are_significant(result):
    up = result.loc[[f"u{i}" for i in range(10)]]
    assert (up["log2FoldChange"] > 0).all()
    assert (up["padj"] < PADJ_MAX).all()


def test_planted_down_genes_negative(result):
    down = result.loc[[f"d{i}" for i in range(10)]]
    assert (down["log2FoldChange"] < 0).all()
    assert (down["padj"] < PADJ_MAX).all()


def test_significant_applies_both_thresholds():
    res = pd.DataFrame({"log2FoldChange": [2.0, 0.5, 2.0, -3.0], "pvalue": 0.0,
                        "padj": [0.01, 0.01, 0.2, 0.001]}, index=list("abcd"))
    assert list(significant(res, PADJ_MAX, LOG2FC_MIN).index) == ["a", "d"]


def test_zero_significant_returns_empty_frame_with_columns():
    res = pd.DataFrame({"log2FoldChange": [0.1], "pvalue": [0.9], "padj": [0.9]}, index=["a"])
    out = significant(res, PADJ_MAX, LOG2FC_MIN)
    assert out.empty and list(out.columns) == list(res.columns)


def test_misaligned_meta_raises_with_missing_ids(planted):
    counts, meta = planted
    with pytest.raises(ValueError, match="P5-cancer"):
        run_dge(counts, meta.drop(index="P5-cancer"))


def test_meta_in_different_order_is_aligned_by_id(planted):
    counts, meta = planted
    shuffled = meta.sample(frac=1, random_state=1)
    a = run_dge(counts, meta)
    b = run_dge(counts, shuffled)
    pd.testing.assert_series_equal(a["log2FoldChange"], b["log2FoldChange"])
