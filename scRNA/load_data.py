import scanpy as sc
import os
from collections import defaultdict
import pandas as pd

def load_data():
    save_path = 'artifacts/scRNA/raw.h5ad'
    data_dir = 'data/scRNA'

    if os.path.exists(save_path):
        print("Loading cached raw data...")
        return sc.read(save_path)
    
    print("No cached data found, building dataset...")

    # Step 1: group files by sample
    samples = defaultdict(dict)

    for f in os.listdir(data_dir):
        if f.endswith('.gz'):
            prefix = f.replace('.matrix.mtx.gz','').replace('.barcodes.tsv.gz','').replace('.genes.tsv.gz','')
            
            if 'matrix' in f:
                samples[prefix]['matrix'] = os.path.join(data_dir, f)
            elif 'barcodes' in f:
                samples[prefix]['barcodes'] = os.path.join(data_dir, f)
            elif 'genes' in f:
                samples[prefix]['genes'] = os.path.join(data_dir, f)

    # Step 2: load each sample
    adatas = []

    for sample, files in samples.items():
        print(f"Loading {sample}")
        
        ad = sc.read_mtx(files['matrix']).T  # transpose is critical
        
        # genes
        genes = pd.read_csv(files['genes'], header=None, sep='\t')
        ad.var_names = genes[1].values
        ad.var_names_make_unique()
        
        # barcodes
        barcodes = pd.read_csv(files['barcodes'], header=None)
        ad.obs_names = barcodes[0].values
        
        # metadata
        ad.obs['sample'] = sample
        
        if 'oLT' in sample:
            ad.obs['condition'] = 'oLT'
        elif 'MT' in sample:
            ad.obs['condition'] = 'MT'
        elif 'SY' in sample:
            ad.obs['condition'] = 'SY'
        
        adatas.append(ad)

    # Step 3: merge
    adata = adatas[0].concatenate(adatas[1:], batch_key='batch')


    # ensure folder exists
    os.makedirs("artifacts/scRNA", exist_ok=True)

    # save
    adata.write(save_path)
    print(f"Saved raw data → {save_path}")

    return adata