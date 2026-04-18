import scanpy as sc
import os

def preprocess(adata):
    save_path = "../artifacts/preprocessed.h5ad"

    # ✅ Load cached version if exists
    if os.path.exists(save_path):
        print("Loading cached preprocessed data...")
        return sc.read(save_path)

    print("Running preprocessing...")

    # Step 1: Ensure unique genes
    adata.var_names_make_unique()

    # Step 2: QC metrics
    print("Computing QC metrics...")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

    # Step 3: Filter cells
    print("Filtering cells...")
    adata = adata[adata.obs.n_genes_by_counts > 500, :]
    adata = adata[adata.obs.n_genes_by_counts < 6000, :]
    adata = adata[adata.obs.pct_counts_mt < 10, :]

    # Step 4: Filter genes
    print("Filtering genes...")
    sc.pp.filter_genes(adata, min_cells=3)

    # Step 5: Normalize
    print("Normalizing...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Step 6: Highly variable genes
    print("Selecting highly variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    adata = adata[:, adata.var.highly_variable]

    # Step 7: Scale + PCA
    print("Scaling and PCA...")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata)

    # save
    os.makedirs("../artifacts", exist_ok=True)
    adata.write(save_path)

    print(f"Saved preprocessed data → {save_path}")
    print(adata)

    return adata