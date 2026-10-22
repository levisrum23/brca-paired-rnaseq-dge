from src import config


def test_thresholds_match_spec():
    assert config.SEED == 42
    assert config.PADJ_MAX == 0.05
    assert config.LOG2FC_MIN == 1.0
    assert config.MAX_DOWNLOAD_MB == 50
    assert config.DATA_DIR.name == "data"
