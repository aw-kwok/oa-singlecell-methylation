import pandas as pd
import os


def map_dmgs():
    dms_path = "../artifacts/methylation/dms.csv"
    annot_path = "../artifacts/methylation/annotation.csv"
    output_path = "../artifacts/methylation/dmgs.csv"

    os.makedirs("../artifacts/methylation", exist_ok=True)

    if os.path.exists(output_path):
        print("dmgs exists, skipping")
        return output_path

    # load data
    dms = pd.read_csv(dms_path)
    annot = pd.read_csv(annot_path)

    # merge CpG → gene
    merged = dms.merge(annot, on="cpg", how="left")

    # drop CpGs without gene mapping
    merged = merged.dropna(subset=["gene_symbol"])

    # aggregate to gene-level
    dmgs = (
        merged
        .groupby("gene_symbol")
        .agg({
            "delta_beta": "mean",   # average methylation change
            "pval": "min"           # strongest signal
        })
        .reset_index()
        .rename(columns={"gene_symbol": "gene"})
    )

    # optional: filter again at gene level
    dmgs = dmgs[
        (dmgs["pval"] < 0.05) &
        (abs(dmgs["delta_beta"]) > 0.1)
    ]

    dmgs.to_csv(output_path, index=False)

    print(f"dmgs created: {len(dmgs)} genes")

    return output_path