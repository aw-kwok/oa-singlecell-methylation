import pandas as pd
import re
import os
from pathlib import Path

def build_beta_matrix():
    input_path = "../data/methylation/GSE73626_non_normalized.txt.gz"
    output_path = "../artifacts/methylation/beta_matrix.csv"

    if os.path.exists(output_path):
        print("beta_matrix exists, skipping")
        return output_path
    
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_path, sep="\t")

    beta = pd.DataFrame(index=df["ID_REF"])

    for col in df.columns:
        if "Methylated Signal" in col:
            sample = col.replace(" Methylated Signal", "")
            u_col = f"{sample} Unmethylated Signal"

            if u_col in df.columns:
                M = df[col]
                U = df[u_col]
                beta[sample] = M / (M + U + 100)

    beta.to_csv(output_path)
    print("beta_matrix created:", beta.shape)

    return output_path