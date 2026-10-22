import pandas as pd

from src.s02_qc import filter_low_counts, pca_separation, run_qc


def test_filter_drops_all_zero_gene(synthetic):
    counts, _ = synthetic
    counts = counts.assign(dead=0)
    assert "dead" not in filter_low_counts(counts).columns


def test_filter_keeps_expressed_gene(synthetic):
    counts, _ = synthetic
    assert "g0" in filter_low_counts(counts).columns


def test_filter_drops_gene_seen_in_too_few_samples(synthetic):
    counts, _ = synthetic
    rare = pd.Series(0, index=counts.index)
    rare.iloc[:2] = 50
    assert "rare" not in filter_low_counts(counts.assign(rare=rare)).columns


def test_pca_separation_positive_for_planted_groups(synthetic):
    counts, meta = synthetic
    assert pca_separation(counts, meta) > 0


def test_run_qc_writes_outputs(synthetic, tmp_path):
    counts, meta = synthetic
    filtered = run_qc(counts, meta, results_dir=tmp_path)
    assert filtered.shape[0] == counts.shape[0]
    for name in ("qc_library_sizes.csv", "qc_pca.png", "qc_correlation.png"):
        assert (tmp_path / name).stat().st_size > 0
    sizes = pd.read_csv(tmp_path / "qc_library_sizes.csv", index_col=0)
    assert {"library_size", "possible_outlier"} <= set(sizes.columns)
