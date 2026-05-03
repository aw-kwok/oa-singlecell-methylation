from methylation.methylation_pipeline import methylation_pipeline
from scRNA.scRNA_pipeline import scRNA_pipeline
from scRNA.integrate import map_dmgs_to_umap

import scanpy as sc


def run_pipeline():
    # 1) run methylation → produces dmgs.csv
    methylation_pipeline()

    # 2) run scRNA → returns adata with UMAP/clusters
    adata = scRNA_pipeline()

    # 3) integrate DMGs onto UMAP
    adata = map_dmgs_to_umap(adata)

    # 4) visualize
    sc.pl.umap(adata, color=["dmg_score", "leiden"])


if __name__ == "__main__":
    run_pipeline()