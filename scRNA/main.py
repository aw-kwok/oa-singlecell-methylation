from load_data import load_data
from preprocess import preprocess
from cluster import cluster

adata = load_data()
adata = preprocess(adata)
adata = cluster(adata)

print(adata)