import gseapy
import pandas as pd

from src.s05_enrichment import run_enrichment, split_up_down


def test_split_up_down_by_sign():
    sig = pd.DataFrame({"log2FoldChange": [2.0, -3.0, 1.5, -1.2]}, index=["a", "b", "c", "d"])
    up, down = split_up_down(sig)
    assert up == ["a", "c"]
    assert down == ["b", "d"]


def test_split_up_down_empty_frame():
    sig = pd.DataFrame({"log2FoldChange": []}, index=[])
    assert split_up_down(sig) == ([], [])


def test_empty_gene_list_returns_empty_frame_without_network(monkeypatch):
    def boom(*a, **k):
        raise AssertionError("network must not be called for an empty gene list")

    monkeypatch.setattr(gseapy, "enrichr", boom)
    out = run_enrichment([], "MSigDB_Hallmark_2020")
    assert out.empty
    assert {"Term", "Adjusted P-value"} <= set(out.columns)


def test_enrichr_uses_library_dict_with_tested_genes_as_background(monkeypatch):
    seen = {}
    library = {"SetA": ["TP53", "MYC"], "SetB": ["EGFR"]}

    class Fake:
        results = pd.DataFrame({"Term": ["t"], "Adjusted P-value": [0.01]})

    monkeypatch.setattr(gseapy, "get_library", lambda name, **k: library if name == "MSigDB_Hallmark_2020" else None)

    def fake_enrichr(**kwargs):
        seen.update(kwargs)
        return Fake()

    monkeypatch.setattr(gseapy, "enrichr", fake_enrichr)
    run_enrichment(["TP53"], "MSigDB_Hallmark_2020", background=["TP53", "MYC", "EGFR", "KRAS"])
    assert seen["gene_sets"] == library
    assert seen["gene_list"] == ["TP53"]
    assert sorted(seen["background"]) == ["EGFR", "KRAS", "MYC", "TP53"]
