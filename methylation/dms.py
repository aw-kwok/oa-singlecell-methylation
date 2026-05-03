import pandas as pd
from scipy.stats import ttest_ind
from pathlib import Path
import numpy as np


def compute_dms():
    beta_path = "artifacts/methylation/beta_matrix.csv"
    meta_path = "artifacts/methylation/sample_metadata.csv"
    output_path = "artifacts/methylation/dms.csv"

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    if Path(output_path).exists():
        print("dms exists, skipping")
        return output_path

    # load data
    beta = pd.read_csv(beta_path, index_col=0)
    meta = pd.read_csv(meta_path)

    # group samples using real metadata
    oa = meta[meta["condition"] == "OA"]["sample_id"].tolist()
    ctrl = meta[meta["condition"] == "control"]["sample_id"].tolist()

    print(f"OA samples: {len(oa)}")
    print(f"Control samples: {len(ctrl)}")

    # sanity check
    if len(oa) < 5 or len(ctrl) < 5:
        raise ValueError("Too few samples in a group")

    beta_oa = beta[oa]
    beta_ctrl = beta[ctrl]

    oa_vals = beta_oa.values
    ctrl_vals = beta_ctrl.values

    # means
    oa_mean = np.nanmean(oa_vals, axis=1)
    ctrl_mean = np.nanmean(ctrl_vals, axis=1)

    delta = oa_mean - ctrl_mean

    # t-test across rows
    stat, pvals = ttest_ind(
        oa_vals,
        ctrl_vals,
        axis=1,
        nan_policy="omit"
    )

    # build dataframe
    dms = pd.DataFrame({
        "cpg": beta.index,
        "delta_beta": delta,
        "pval": pvals
    })

    # filter strong signals
    dms = dms[
        (np.abs(dms["delta_beta"]) > 0.2) &
        (dms["pval"] < 0.05)
    ]

    # direction
    dms["direction"] = np.where(
        dms["delta_beta"] < 0,
        "hypo",
        "hyper"
    )

    dms.to_csv(output_path, index=False)

    print(f"dms created: {len(dms)} CpGs")

    return output_path