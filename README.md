# OA Multi-Modal Analysis (scRNA-seq + DNA Methylation)

## Overview

This repository provides a reproducible pipeline to analyze osteoarthritis (OA) using:

1. **Single-cell RNA-seq (scRNA-seq)** to resolve disease-associated vs homeostatic cell states (UMAP + Leiden clustering)
2. **Bulk DNA methylation (Illumina 450k)** to identify differentially methylated sites (DMS) and genes (DMGs)
3. **Integration** by projecting methylation-derived gene sets onto the scRNA embedding via a per-cell gene score (`dmg_score`)

The code is organized as lightweight, cache-aware scripts (write to `artifacts/`) plus notebooks for interpretation and score construction.

---

## Data

- **scRNA-seq**: GEO **GSE152805** (cartilage and synovium single-cell expression)
- **DNA methylation**: GEO **GSE73626** (450k array; used by the methylation pipeline scripts)

---

## Project Structure

```text
oa-singlecell-methylation/

├── main.py                        # end-to-end run: methylation -> scRNA -> DMG overlay
├── requirements.txt
├── README.md

├── data/
│   ├── scRNA/                     # raw per-sample matrices (.gz)
│   │   ├── GSM4626763_SY_113.matrix.mtx.gz
│   │   ├── GSM4626763_SY_113.barcodes.tsv.gz
│   │   ├── GSM4626763_SY_113.genes.tsv.gz
│   │   ├── GSM4626764_SY_116.matrix.mtx.gz
│   │   ├── GSM4626764_SY_116.barcodes.tsv.gz
│   │   ├── GSM4626764_SY_116.genes.tsv.gz
│   │   ├── ...
│   │   ├── GSM4626771_OA_MT_118.matrix.mtx.gz
│   │   ├── GSM4626771_OA_MT_118.barcodes.tsv.gz
│   │   └── GSM4626771_OA_MT_118.genes.tsv.gz
│   └── methylation/               # raw GSE73626 + 450k manifest
│       ├── GSE73626_non_normalized.txt.gz
│       ├── GSE73626_series_matrix.txt.gz
│       └── humanmethylation450_15017482_v1-2.csv

├── scRNA/
│   ├── load_data.py               # read .mtx/.genes/.barcodes per sample, concatenate
│   ├── preprocess.py              # QC, filtering, normalization, HVGs, PCA
│   ├── cluster.py                 # neighbors, UMAP, Leiden
│   ├── integrate.py               # DMG gene-set scoring on UMAP (writes clustered_with_dmgs.h5ad)
│   └── scRNA_pipeline.py          # script entrypoint for scRNA pipeline

├── methylation/
│   ├── build_beta.py              # build beta matrix + sample metadata from GSE73626
│   ├── annotation.py              # parse 450k manifest -> CpG -> gene mapping
│   ├── dms.py                     # t-test per CpG to call DMS
│   ├── dmgs.py                    # aggregate DMS -> DMGs
│   └── methylation_pipeline.py    # script entrypoint for methylation pipeline

├── notebooks/
│   ├── scRNA_analysis.ipynb        # interpretation + OA score (writes artifacts/scRNA/oa_scores.csv)
│   └── methylation_overlay.ipynb   # visualize methylation overlay on scRNA UMAP

└── artifacts/                      # generated outputs (safe to delete/rebuild)
    ├── scRNA/
    └── methylation/
```

---

## Setup

```bash
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt
```

---

## Required Input Files

### scRNA-seq inputs

Place the (gzipped) MatrixMarket + manifest files under `data/scRNA/`.
Each sample must have three files using the exact suffixes shown below:

```text
data/scRNA/
  GSM4626763_SY_113.matrix.mtx.gz
  GSM4626763_SY_113.barcodes.tsv.gz
  GSM4626763_SY_113.genes.tsv.gz
  ...
```

The loader groups files by the shared sample prefix (e.g., `GSM4626763_SY_113`) and builds one AnnData per sample before concatenation.

### Methylation inputs

The methylation scripts expect the following files:

```text
data/methylation/
  GSE73626_non_normalized.txt.gz
  GSE73626_series_matrix.txt.gz
  humanmethylation450_15017482_v1-2.csv
```

---

## Running the Pipelines

### Option A: Run the full end-to-end pipeline

From the project root:

```bash
python main.py
```

This will:

1. Run methylation processing and write `artifacts/methylation/dmgs.csv`
2. Run scRNA preprocessing + clustering and write `artifacts/scRNA/clustered.h5ad`
3. Compute `dmg_score` on the scRNA object and write `artifacts/scRNA/clustered_with_dmgs.h5ad`

### Option B: Run scRNA-only

```bash
python scRNA/scRNA_pipeline.py
```

### Option C: Run methylation-only

```bash
python methylation/methylation_pipeline.py
```

---

## Notebooks

- `notebooks/scRNA_analysis.ipynb`
  - visualizes clusters and condition shifts
  - constructs a sample-level OA score
  - writes: `artifacts/scRNA/oa_scores.csv`

- `notebooks/methylation_overlay.ipynb`
  - loads clustered scRNA data and methylation DMGs
  - visualizes DMG-derived scores (e.g., `dmg_score`) on the UMAP embedding

---

## Outputs (Artifacts)

### scRNA

- `artifacts/scRNA/raw.h5ad`
- `artifacts/scRNA/preprocessed.h5ad`
- `artifacts/scRNA/clustered.h5ad`
- `artifacts/scRNA/clustered_with_dmgs.h5ad`
- `artifacts/scRNA/oa_scores.csv` (created by the scRNA notebook)

### Methylation

- `artifacts/methylation/beta_matrix.csv`
- `artifacts/methylation/sample_metadata.csv`
- `artifacts/methylation/annotation.csv`
- `artifacts/methylation/dms.csv`
- `artifacts/methylation/dmgs.csv`

---

## Notes

- Most steps are cache-aware: if the expected output file exists in `artifacts/`, the step will print a message and skip recomputation.
- If you change inputs and want a clean rebuild, delete the relevant outputs under `artifacts/`.
