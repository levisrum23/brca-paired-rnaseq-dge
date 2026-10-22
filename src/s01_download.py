import hashlib
import urllib.request
from pathlib import Path


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
