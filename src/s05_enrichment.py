import gseapy as gp
import pandas as pd

from src.config import PADJ_MAX, RESULTS_DIR

LIBRARIES = {"hallmark": "MSigDB_Hallmark_2020", "kegg": "KEGG_2026"}
RESULT_COLUMNS = ["Gene_set", "Term", "Overlap", "P-value", "Adjusted P-value", "Odds Ratio",
                  "Combined Score", "Genes"]


def split_up_down(sig: pd.DataFrame) -> tuple[list[str], list[str]]:
    """Gene lists higher (log2FC > 0) and lower (log2FC < 0) in tumor."""
    return (list(sig.index[sig["log2FoldChange"] > 0]), list(sig.index[sig["log2FoldChange"] < 0]))


def run_enrichment(genes: list[str], library: str, background: list[str] | None = None) -> pd.DataFrame:
    """Over-representation test (human). Empty gene list gives an empty frame, no network call.

    With a background (the genes that were tested), the library is fetched once and the test runs
    locally against that background. Without one, Enrichr's default background is used.
    """
    if not genes:
        return pd.DataFrame(columns=RESULT_COLUMNS)
    if background is None:
        return gp.enrichr(gene_list=genes, gene_sets=library, organism="human", outdir=None).results
    gene_sets = gp.get_library(library, organism="Human")
    return gp.enrichr(gene_list=genes, gene_sets=gene_sets, background=list(background),
                      organism="human", outdir=None).results


def main() -> None:
    sig = pd.read_csv(RESULTS_DIR / "dge_significant.csv", index_col=0)
    background = list(pd.read_csv(RESULTS_DIR / "dge_all.csv", index_col=0).index)
    up, down = split_up_down(sig)
    for direction, genes in (("up", up), ("down", down)):
        for key, library in LIBRARIES.items():
            res = run_enrichment(genes, library, background)
            res.to_csv(RESULTS_DIR / f"enrichment_{direction}_{key}.csv", index=False)
            n_sig = int((res["Adjusted P-value"] < PADJ_MAX).sum()) if len(res) else 0
            print(f"{direction}/{key}: {len(genes)} genes, {n_sig} terms with adj. p < 0.05")


if __name__ == "__main__":
    main()
