import hashlib
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from src.config import DATA_DIR, MAX_DOWNLOAD_MB

ANNOTATION_COLS = ["Chr", "Start", "End", "Strand", "Length"]

# GEO series GSE133998. The checksum was recorded on first download (GEO publishes none).
COUNTS_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133998/suppl/GSE133998_cancer_and_adjacent_normal.txt.gz"
COUNTS_SHA256 = "3983a2180e3631bfc7141c7fd0241e9ae12679444b4a7012df583ec4ec08f590"
COUNTS_PATH = DATA_DIR / "raw" / "GSE133998_cancer_and_adjacent_normal.txt.gz"
META_PATH = DATA_DIR / "metadata.csv"


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download_verified(url: str, dest: Path, expected_sha256: str | None, max_mb: int) -> str:
    """Download url to dest. Delete dest and raise RuntimeError on size or checksum failure."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    limit = max_mb * 1024 * 1024
    written = 0
    try:
        with urllib.request.urlopen(url) as resp, open(dest, "wb") as out:
            for chunk in iter(lambda: resp.read(1 << 20), b""):
                written += len(chunk)
                if written > limit:
                    raise RuntimeError(f"download exceeds {max_mb} MB limit")
                out.write(chunk)
        actual = sha256_of(dest)
        if expected_sha256 is not None and actual != expected_sha256:
            raise RuntimeError(f"checksum mismatch: expected {expected_sha256}, got {actual}")
        return actual
    except BaseException:
        dest.unlink(missing_ok=True)
        raise


def load_dataset(counts_path: Path, meta_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (counts, meta): counts is samples x genes (int); meta is indexed by sample ID."""
    raw = pd.read_csv(counts_path, sep="\t", index_col=0)
    raw = raw.drop(columns=[c for c in ANNOTATION_COLS if c in raw.columns])
    blank = raw.isna().all(axis=1)
    if blank.any():
        print(f"dropped {int(blank.sum())} gene rows with no values: {list(raw.index[blank])}")
        raw = raw[~blank]
    if raw.isna().any().any():
        raise ValueError("counts contain missing values in partly filled rows")
    values = raw.to_numpy(dtype=float)
    if not np.array_equal(values, np.round(values)):
        raise ValueError("counts are not integers (TPM/FPKM/normalized data cannot be used by DESeq2)")
    n_dupes = int(raw.index.duplicated().sum())
    if n_dupes:
        print(f"merged {n_dupes} duplicate gene rows by summing")
    counts = raw.astype("int64").groupby(level=0, sort=False).sum().T

    meta = pd.read_csv(meta_path, index_col="sample")
    missing = [s for s in counts.index if s not in meta.index]
    if missing:
        raise ValueError(f"samples missing from metadata: {missing}")
    return counts, meta.loc[counts.index]


def main() -> None:
    if not COUNTS_PATH.exists():
        download_verified(COUNTS_URL, COUNTS_PATH, COUNTS_SHA256, MAX_DOWNLOAD_MB)
    elif sha256_of(COUNTS_PATH) != COUNTS_SHA256:
        raise RuntimeError(f"{COUNTS_PATH} exists but its checksum does not match")
    counts, meta = load_dataset(COUNTS_PATH, META_PATH)
    print(f"ok: {counts.shape[0]} samples x {counts.shape[1]} genes; {meta['condition'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
