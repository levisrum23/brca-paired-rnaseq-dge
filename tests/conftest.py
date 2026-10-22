import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def synthetic():
    """12 samples (6 tumor, 6 normal, paired) x 30 genes; genes g0-g4 are high in tumor."""
    rng = np.random.default_rng(42)
    samples = [f"P{i}-{t}" for i in range(6) for t in ("adj", "cancer")]
    base = rng.negative_binomial(5, 0.02, size=(12, 30))
    tumor = np.array([s.endswith("cancer") for s in samples])
    base[np.ix_(tumor, range(5))] *= 6
    counts = pd.DataFrame(base, index=samples, columns=[f"g{i}" for i in range(30)])
    meta = pd.DataFrame(
        {"condition": np.where(tumor, "tumor", "normal"), "patient": [s.split("-")[0] for s in samples]},
        index=samples,
    )
    return counts, meta
