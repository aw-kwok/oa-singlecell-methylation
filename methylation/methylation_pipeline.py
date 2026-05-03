from build_beta import build_beta_matrix
from annotation import build_annotation
from dms import compute_dms
from dmgs import map_dmgs


def methylation_pipeline():
    beta = build_beta_matrix()
    annotation = build_annotation()
    dms = compute_dms()
    dmgs = map_dmgs()

if __name__ == "__main__":
    methylation_pipeline()