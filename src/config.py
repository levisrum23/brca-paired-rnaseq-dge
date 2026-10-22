from pathlib import Path

SEED = 42
PADJ_MAX = 0.05
LOG2FC_MIN = 1.0
MAX_DOWNLOAD_MB = 50

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
