import gseapy
import numpy as np
import pandas as pd
import pytest

from src.config import DATA_DIR, SEED
from run_all import check_sanity, run_pipeline

EXPECTED = [
    "qc_library_sizes.csv", "qc_summary.csv", "qc_pca.png", "qc_correlation.png",
    "dge_all.csv", "dge_significant.csv",
    "volcano.png", "ma_plot.png", "top_genes_heatmap.png",
    "enrichment_up_hallmark.csv", "enrichment_up_kegg.csv",
    "enrichment_down_hallmark.csv", "enrichment_down_kegg.csv",
    "marker_comparison.csv",
]


@pytest.fixture
def dataset(tmp_path):
    """Small GEO-style counts file plus metadata: 6 pairs x 60 genes, 10 planted up and 10 down."""
    rng = np.random.default_rng(SEED)
    samples = [f"P{i}-{t}" for i in range(6) for t in ("adj", "cancer")]
    tumor = np.array([s.endswith("cancer") for s in samples])
    mu = np.full((12, 60), 500.0)
    mu[np.ix_(tumor, range(10))] *= 8
    mu[np.ix_(tumor, range(10, 20))] /= 8
    genes = [f"u{i}" for i in range(10)] + [f"d{i}" for i in range(10)] + [f"n{i}" for i in range(40)]
    table = pd.DataFrame(rng.poisson(mu).T, index=genes, columns=samples)
    table.insert(0, "Chr", "chr1")
    table.insert(1, "Start", 1)
    table.insert(2, "End", 2)
    table.insert(3, "Strand", "+")
    table.insert(4, "Length", 100)
    table.index.name = "Geneid"
    counts_path = tmp_path / "counts.tsv"
    table.to_csv(counts_path, sep="\t")
    meta = pd.DataFrame({"sample": samples, "condition": np.where(tumor, "tumor", "normal"),
                         "patient": [s.split("-")[0] for s in samples]})
    meta_path = tmp_path / "meta.csv"
    meta.to_csv(meta_path, index=False)
    return counts_path, meta_path


def test_rerun_from_clean_results(dataset, tmp_path, monkeypatch):
    class Fake:
        results = pd.DataFrame({"Gene_set": ["x"], "Term": ["t"], "Overlap": ["1/2"], "P-value": [0.01],
                                "Adjusted P-value": [0.02], "Odds Ratio": [1.0], "Combined Score": [1.0],
                                "Genes": ["u0"]})

    monkeypatch.setattr(gseapy, "get_library", lambda name, **k: {"SetA": ["u0", "u1"]})
    monkeypatch.setattr(gseapy, "enrichr", lambda **k: Fake())
    results = tmp_path / "fresh_results"
    assert not results.exists()
    run_pipeline(*dataset, markers_path=DATA_DIR / "known_markers.csv", results_dir=results)
    missing = [f for f in EXPECTED if not (results / f).exists() or (results / f).stat().st_size == 0]
    assert missing == []


def test_sanity_check_fails_when_groups_do_not_separate():
    with pytest.raises(RuntimeError, match="PCA"):
        check_sanity(pca_sep=-0.1, n_significant=10)


def test_sanity_check_warns_but_passes_when_nothing_significant(capsys):
    check_sanity(pca_sep=0.3, n_significant=0)
    assert "warning" in capsys.readouterr().out.lower()
