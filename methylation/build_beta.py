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
    df = df.set_index("ID_REF")

    beta = pd.DataFrame(index=df.index)

    for col in df.columns:
        if "Unmethylated Signal" in col:
            sample = col.replace(" Unmethylated Signal", "")
            m_col = f"{sample} Methylated Signal"

            if m_col in df.columns:
                U = df[col]
                M = df[m_col]

                beta[sample] = M / (M + U + 100)
            else:
                print(f"Skipping {sample} (missing Methylated column)")
    
    beta.to_csv(output_path)
    print("beta_matrix created:", beta.shape)

    return output_path