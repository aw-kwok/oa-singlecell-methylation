from methylation.build_beta import build_beta_matrix
from methylation.annotation import build_annotation
from methylation.dms import compute_dms
from methylation.dmgs import map_dmgs


def methylation_pipeline():
    beta = build_beta_matrix()
    annotation = build_annotation()
    dms = compute_dms()
    dmgs = map_dmgs()

if __name__ == "__main__":
    methylation_pipeline()