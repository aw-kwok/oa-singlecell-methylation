import scanpy as sc
import os

def cluster(adata):
    save_path = "artifacts/scRNA/clustered.h5ad"

    # Step 0: Load cached clustered data if it exists
    if os.path.exists(save_path):
        print("Loading cached clustered data...")
        return sc.read(save_path)

    print("Running clustering...")

    # Step 1: Compute neighborhood graph using PCA representation
    print("Computing neighbors...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

    # Step 2: Compute UMAP embedding for visualization
    print("Computing UMAP...")
    sc.tl.umap(adata)

    # Step 3: Perform Leiden clustering to identify cell states
    print("Running Leiden clustering...")
    sc.tl.leiden(adata, resolution=0.5)

    # Step 4: Save clustered AnnData object for reuse
    os.makedirs("artifacts/scRNA", exist_ok=True)
    adata.write(save_path)

    print(f"Saved clustered data → {save_path}")
    print(adata)

    return adata