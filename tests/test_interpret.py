import pandas as pd

from src.s06_interpret import compare_markers


def _res():
    return pd.DataFrame({"log2FoldChange": [2.5, -1.5, 0.3, float("nan")],
                         "padj": [0.001, 0.02, 0.9, float("nan")]}, index=["UP1", "DOWN1", "FLAT", "NANGENE"])


def _markers(rows):
    return pd.DataFrame(rows, columns=["gene", "expected_direction", "source"])


def test_agrees_true_when_direction_matches():
    out = compare_markers(_res(), _markers([("UP1", "up", "s"), ("DOWN1", "down", "s")]))
    assert list(out["agrees"]) == [True, True]
    assert list(out["observed_lfc"]) == [2.5, -1.5]


def test_agrees_false_when_opposite():
    out = compare_markers(_res(), _markers([("UP1", "down", "s")]))
    assert list(out["agrees"]) == [False]


def test_missing_gene_gives_none():
    out = compare_markers(_res(), _markers([("ABSENT", "up", "s")]))
    assert out.loc[0, "agrees"] is None


def test_gene_without_result_gives_none():
    out = compare_markers(_res(), _markers([("NANGENE", "up", "s")]))
    assert out.loc[0, "agrees"] is None


def test_non_significant_marker_keeps_padj_visible():
    out = compare_markers(_res(), _markers([("FLAT", "up", "s")]))
    assert out.loc[0, "padj"] == 0.9
