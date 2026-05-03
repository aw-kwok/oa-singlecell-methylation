import pandas as pd
import scanpy as sc


def map_dmgs_to_umap(adata):
    # load DMGs
    dmgs = pd.read_csv("artifacts/methylation/dmgs.csv")
    genes = dmgs["gene"].tolist()

    # match genes to scRNA
    genes = [g for g in genes if g in adata.var_names]

    print(f"Matched genes: {len(genes)}")

    if len(genes) == 0:
        raise ValueError("No DMGs found in scRNA dataset")

    # compute gene score
    sc.tl.score_genes(
        adata,
        gene_list=genes,
        score_name="dmg_score"
    )

    adata.write("artifacts/scRNA/clustered_with_dmgs.h5ad")

    return adata