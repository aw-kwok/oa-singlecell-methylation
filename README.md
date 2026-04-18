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

## Requirements

This project uses Python and the following core libraries:

- `scanpy`
- `anndata`
- `pandas`
- `numpy`
- `scipy`
- `matplotlib`
- `seaborn`

To install dependencies:

```bash
pip install -r requirements.txt