import pandas as pd
from pathlib import Path
import gzip
import re


def parse_sample_titles(matrix_path):
    if str(matrix_path).endswith(".gz"):
        f = gzip.open(matrix_path, mode="rt", encoding="utf-8", errors="ignore")
    else:
        f = open(matrix_path, mode="r", encoding="utf-8", errors="ignore")

    with f:
        for line in f:
            if line.startswith("!Sample_title"):
                return re.findall(r'"(.*?)"', line)

    raise ValueError("!Sample_title not found")


def build_beta_matrix():
    input_path="data/methylation/GSE73626_non_normalized.txt.gz"
    matrix_path="data/methylation/GSE73626_series_matrix.txt.gz"
    beta_output="artifacts/methylation/beta_matrix.csv"
    metadata_output="artifacts/methylation/sample_metadata.csv"
    
    Path(beta_output).parent.mkdir(parents=True, exist_ok=True)
    Path(metadata_output).parent.mkdir(parents=True, exist_ok=True)

    if Path(beta_output).exists() and Path(metadata_output).exists():
        print("beta + metadata exist, skipping")
        return beta_output

    # load data
    df = pd.read_csv(input_path, sep="\t")
    df = df.set_index("ID_REF")

    # parse GEO metadata
    titles = parse_sample_titles(matrix_path)

    beta = pd.DataFrame(index=df.index)

    metadata = []
    dropped = []

    # build beta + track samples
    for col in df.columns:
        if "Unmethylated Signal" in col:
            sample = col.replace(" Unmethylated Signal", "")
            m_col = f"{sample} Methylated Signal"

            if m_col in df.columns:
                U = df[col]
                M = df[m_col]

                beta[sample] = M / (M + U + 100)
                metadata.append(sample)
            else:
                print(f"Skipping {sample} (missing M column)")
                dropped.append(sample)

    print(f"Beta shape: {beta.shape}")
    print(f"Dropped samples: {dropped}")

    # align titles (safe because order preserved)
    titles = titles[:len(metadata)]

    meta_df = pd.DataFrame({
        "sample_id": metadata,
        "title": titles
    })

    # derive labels
    meta_df["condition"] = meta_df["title"].apply(
        lambda x: "control" if "control" in x.lower() else "OA"
    )

    meta_df["joint"] = meta_df["title"].apply(
        lambda x: "hip" if "hip" in x.lower()
        else ("knee" if "knee" in x.lower() else None)
    )

    # save outputs
    beta.to_csv(beta_output)
    meta_df.to_csv(metadata_output, index=False)

    print("beta_matrix + metadata created")

    return beta_output