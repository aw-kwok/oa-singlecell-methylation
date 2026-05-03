import pandas as pd
import os
from pathlib import Path

def build_annotation():
    input_path = "../data/methylation/humanmethylation450_15017482_v1-2.csv"
    output_path = "../artifacts/methylation/annotation.csv"

    if os.path.exists(output_path):
        print("annotation exists, skipping")
        return output_path
    
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(
        input_path,
        low_memory=False,
        skiprows=7
    )

    annot = df[["Name", "UCSC_RefGene_Name"]].copy()
    annot.columns = ["cpg", "gene_symbol"]

    annot["gene_symbol"] = annot["gene_symbol"].str.split(";").str[0]
    annot = annot.dropna()

    annot.to_csv(output_path, index=False)

    print("annotation created:", annot.shape)
    return output_path