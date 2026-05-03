from build_beta import build_beta_matrix
from annotation import build_annotation
from dms import compute_dms

beta = build_beta_matrix()
annotation = build_annotation()

compute_dms()