import gseapy as gp
import pandas as pd

from src.config import RESULTS_DIR

LIBRARIES = {"hallmark": "MSigDB_Hallmark_2020", "kegg": "KEGG_2026"}
RESULT_COLUMNS = ["Gene_set", "Term", "Overlap", "P-value", "Adjusted P-value", "Odds Ratio",
                  "Combined Score", "Genes"]


def split_up_down(sig: pd.DataFrame) -> tuple[list[str], list[str]]:
    """Gene lists higher (log2FC > 0) and lower (log2FC < 0) in tumor."""
    return (list(sig.index[sig["log2FoldChange"] > 0]), list(sig.index[sig["log2FoldChange"] < 0]))


def run_enrichment(genes: list[str], library: str) -> pd.DataFrame:
    """Over-representation test via Enrichr (human). Empty gene list gives an empty frame, no network call."""
    if not genes:
        return pd.DataFrame(columns=RESULT_COLUMNS)
    return gp.enrichr(gene_list=genes, gene_sets=library, organism="human", outdir=None).results


def main() -> None:
    sig = pd.read_csv(RESULTS_DIR / "dge_significant.csv", index_col=0)
    up, down = split_up_down(sig)
    for direction, genes in (("up", up), ("down", down)):
        for key, library in LIBRARIES.items():
            res = run_enrichment(genes, library)
            res.to_csv(RESULTS_DIR / f"enrichment_{direction}_{key}.csv", index=False)
            n_sig = int((res["Adjusted P-value"] < 0.05).sum()) if len(res) else 0
            print(f"{direction}/{key}: {len(genes)} genes, {n_sig} terms with adj. p < 0.05")


if __name__ == "__main__":
    main()
