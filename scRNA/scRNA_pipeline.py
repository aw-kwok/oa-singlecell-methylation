from load_data import load_data
from preprocess import preprocess
from cluster import cluster

def scRNA_pipeline():
    adata = load_data()
    adata = preprocess(adata)
    adata = cluster(adata)

    print(adata)

if __name__ == "__main__":
    scRNA_pipeline()