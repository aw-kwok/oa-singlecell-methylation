from load_data import load_data
from preprocess import preprocess

adata = load_data()
adata = preprocess(adata)

print(adata)