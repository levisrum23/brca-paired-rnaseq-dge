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
