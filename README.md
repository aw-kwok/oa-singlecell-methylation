# OA Multi-Modal Analysis

## Overview

This project analyzes human osteoarthritis (OA) using multiple genomic data modalities, beginning with single-cell RNA sequencing (scRNA-seq) and extending to DNA methylation.

The central goal is to understand how osteoarthritis alters cellular composition and regulatory programs within joint tissues. Rather than treating OA as a binary condition, this project models it as a shift in underlying cell states.

The workflow proceeds in two stages:

1. **scRNA-seq analysis** to identify distinct cell states and characterize how their proportions change across tissue conditions
2. **Methylation analysis (planned)** to determine whether epigenetic signals capture the same disease-associated shifts

By integrating these modalities, the project aims to link gene expression changes with upstream regulatory mechanisms.

---

## Data

- **scRNA-seq**: GSE152805  
  Human osteoarthritic cartilage and synovium at single-cell resolution

- **DNA methylation (planned)**: GEO datasets such as GSE63695  

---

## Project Structure

The repository is organized into raw data, processing pipelines, and analysis outputs:

```text
oa-singlecell-methylation/

├── data/                          # raw GEO scRNA-seq data
│   └── scRNA/
│       ├── GSM4626763_SY_113.*
│       ├── GSM4626764_SY_116.*
│       ├── GSM4626765_SY_118.*
│       ├── GSM4626766_OA_oLT_113.*
│       ├── GSM4626767_OA_oLT_116.*
│       ├── GSM4626768_OA_oLT_118.*
│       ├── GSM4626769_OA_MT_113.*
│       ├── GSM4626770_OA_MT_116.*
│       └── GSM4626771_OA_MT_118.*

├── artifacts/                     # processed outputs
│   └── scRNA/
│       ├── raw.h5ad               # merged raw dataset
│       ├── preprocessed.h5ad      # QC + normalized data
│       ├── clustered.h5ad         # PCA, UMAP, Leiden clusters
│       └── oa_scores.csv          # sample-level OA scores

├── scRNA/                         # scRNA processing pipeline
│   ├── load_data.py               # load and merge GEO matrices
│   ├── preprocess.py              # QC, filtering, normalization, PCA
│   ├── cluster.py                 # neighbors, UMAP, Leiden clustering
│   └── main.py                    # end-to-end pipeline runner

├── notebooks/                     # analysis and interpretation
│   └── scRNA_analysis.ipynb       # UMAP, marker genes, OA score

├── requirements.txt               # Python dependencies
└── README.md                      # project documentation
```

---

## Requirements

This project uses Python and the following core libraries:

- `scanpy`
- `anndata`
- `pandas`
- `numpy`
- `scipy`
- `matplotlib`
- `seaborn`
- `igraph`
- `leidenalg`
- `ipykernel`

To install dependencies:

```bash
pip install -r requirements.txt
```

---

## How to Run

## How to Run

### 1. Setup environment

Create and activate a virtual environment, then install dependencies:

```bash
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt
```

---

### 2. Prepare data

Place the raw scRNA-seq files (GEO download) in:

```text
data/scRNA/
```

The pipeline expects `.mtx`, `.genes.tsv.gz`, and `.barcodes.tsv.gz` files for each sample.

---

### 3. Run the pipeline

Execute the full scRNA processing pipeline:

```bash
python scRNA/main.py
```

---

### 4. What the pipeline does

The pipeline will:

1. Load raw data (`load_data.py`)
2. Preprocess and normalize (`preprocess.py`)
3. Perform clustering (PCA, UMAP, Leiden) (`cluster.py`)
4. Save outputs to:

```text
artifacts/scRNA/
```

including:
- `raw.h5ad`
- `preprocessed.h5ad`
- `clustered.h5ad`

---

### 5. Run analysis notebook

Open and run:

```text
notebooks/scRNA_analysis.ipynb
```

This notebook:
- visualizes UMAP clusters
- performs marker gene analysis
- computes the OA score

---

### 6. Output

Final results include:

```text
artifacts/scRNA/oa_scores.csv
```

which contains sample-level OA scores used for downstream analysis.

---

### Notes

- If processed `.h5ad` files already exist, the pipeline will load cached data instead of recomputing.
- Ensure you are running commands from the project root directory.

---

## To Do

Future work will extend this analysis by integrating DNA methylation data to predict the RNA-derived OA score.

This will allow us to test whether epigenetic variation can capture and predict disease-associated shifts in cellular composition.