import pytest

from src.s01_download import download_verified, sha256_of

ABC_SHA256 = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"


def _source(tmp_path, payload=b"abc"):
    src = tmp_path / "source.bin"
    src.write_bytes(payload)
    return src.as_uri()


def test_sha256_of_known_bytes(tmp_path):
    f = tmp_path / "x.bin"
    f.write_bytes(b"abc")
    assert sha256_of(f) == ABC_SHA256


def test_download_returns_actual_sha256_on_match(tmp_path):
    dest = tmp_path / "out.bin"
    assert download_verified(_source(tmp_path), dest, ABC_SHA256, max_mb=1) == ABC_SHA256
    assert dest.read_bytes() == b"abc"


def test_download_removes_file_on_checksum_mismatch(tmp_path):
    dest = tmp_path / "out.bin"
    with pytest.raises(RuntimeError):
        download_verified(_source(tmp_path), dest, "0" * 64, max_mb=1)
    assert not dest.exists()


def test_download_rejects_over_size(tmp_path):
    dest = tmp_path / "out.bin"
    with pytest.raises(RuntimeError):
        download_verified(_source(tmp_path), dest, None, max_mb=0)
    assert not dest.exists()


def _write_counts(tmp_path, rows, samples=("P1-adj", "P1-cancer")):
    p = tmp_path / "counts.tsv"
    header = "Geneid\tChr\tStart\tEnd\tStrand\tLength\t" + "\t".join(samples)
    lines = [header] + [f"{g}\tchr1\t1\t2\t+\t100\t" + "\t".join(str(v) for v in vals) for g, vals in rows]
    p.write_text("\n".join(lines) + "\n")
    return p


def _write_meta(tmp_path, ids=("P1-adj", "P1-cancer")):
    p = tmp_path / "meta.csv"
    cond = {"adj": "normal", "cancer": "tumor"}
    lines = ["sample,condition,patient"] + [f"{i},{cond[i.split('-')[1]]},{i.split('-')[0]}" for i in ids]
    p.write_text("\n".join(lines) + "\n")
    return p


def test_load_returns_samples_by_genes_with_condition(tmp_path):
    from src.s01_download import load_dataset
    counts, meta = load_dataset(_write_counts(tmp_path, [("A", (1, 5)), ("B", (2, 8))]), _write_meta(tmp_path))
    assert list(counts.index) == ["P1-adj", "P1-cancer"]
    assert list(counts.columns) == ["A", "B"]
    assert list(meta["condition"]) == ["normal", "tumor"]


def test_load_rejects_non_integer_counts(tmp_path):
    from src.s01_download import load_dataset
    with pytest.raises(ValueError, match="counts are not integers"):
        load_dataset(_write_counts(tmp_path, [("A", (1.5, 5))]), _write_meta(tmp_path))


def test_load_merges_duplicate_genes(tmp_path):
    from src.s01_download import load_dataset
    counts, _ = load_dataset(_write_counts(tmp_path, [("A", (1, 5)), ("A", (2, 3))]), _write_meta(tmp_path))
    assert list(counts.columns) == ["A"]
    assert list(counts["A"]) == [3, 8]


def test_load_errors_on_missing_metadata_ids(tmp_path):
    from src.s01_download import load_dataset
    with pytest.raises(ValueError, match="P1-cancer"):
        load_dataset(_write_counts(tmp_path, [("A", (1, 5))]), _write_meta(tmp_path, ids=("P1-adj",)))


def test_load_drops_fully_blank_gene_rows(tmp_path):
    from src.s01_download import load_dataset
    p = _write_counts(tmp_path, [("A", (1, 5)), ("B", (2, 8))])
    p.write_text(p.read_text() + "BLANK\tchr1\t1\t2\t+\t100\t\t\n")
    counts, _ = load_dataset(p, _write_meta(tmp_path))
    assert list(counts.columns) == ["A", "B"]


def test_load_rejects_partly_blank_gene_rows(tmp_path):
    from src.s01_download import load_dataset
    p = _write_counts(tmp_path, [("A", (1, 5))])
    p.write_text(p.read_text() + "PART\tchr1\t1\t2\t+\t100\t3\t\n")
    with pytest.raises(ValueError, match="missing values"):
        load_dataset(p, _write_meta(tmp_path))
