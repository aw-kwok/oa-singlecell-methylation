from scRNA.load_data import load_data
from scRNA.preprocess import preprocess
from scRNA.cluster import cluster

def scRNA_pipeline():
    adata = load_data()
    adata = preprocess(adata)
    adata = cluster(adata)

    print(adata)

    return adata

if __name__ == "__main__":
    scRNA_pipeline()